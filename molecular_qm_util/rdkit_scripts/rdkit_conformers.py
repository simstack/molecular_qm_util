# conformer_generator.py
from __future__ import annotations

import concurrent.futures
import logging
import os
import time
from pathlib import Path

from simstack.core.node import node
from molecular_qm_models import Atom, Molecule, MoleculeList
from molecular_qm_models.prune_conformers import prune_conformers
from simstack.models import StringData

from molecular_qm_util.obabel_scripts.openbabel_conformers import ConformerGenerationInput

logger = logging.getLogger(__name__)

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem
except Exception as e:  # pragma: no cover
    Chem = None  # type: ignore[assignment]
    AllChem = None  # type: ignore[assignment]
    _RDKIT_IMPORT_ERROR = e
else:
    _RDKIT_IMPORT_ERROR = None

# Default RMSD threshold (Angstrom) used when pruning RDKit conformers.
_DEFAULT_PRUNE_RMS_THRESH = 0.1


def _require_rdkit() -> None:
    if Chem is None or AllChem is None:
        raise RuntimeError(
            "RDKit is not available in this environment. Install it first (e.g. `pip install rdkit-pypi`) "
            "or use --method crest with an external CREST installation."
        ) from _RDKIT_IMPORT_ERROR


def _rdkit_conf_to_molecule(mol, conf, energy: float, smiles) -> Molecule:
    """Convert a single RDKit conformer into a qm_models Molecule."""
    new_molecule = Molecule()
    new_molecule.smiles = smiles
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        new_molecule.add_atom(
            Atom(
                element=atom.GetSymbol(),
                x=pos.x,
                y=pos.y,
                z=pos.z,
            )
        )
    new_molecule.properties["energy"] = energy
    return new_molecule


@node
def conformers_rdkit(molecule: Molecule, confgen_input: ConformerGenerationInput, **kwargs) -> MoleculeList:
    node_runner = kwargs["node_runner"]
    num_confs = confgen_input.num_confs
    seed = confgen_input.seed
    profile = confgen_input.profile
    molecule_smiles = molecule.smiles if molecule.smiles is not None else "NA"
    molecule_formula = molecule.formula if molecule.formula is not None else "NA"
    node_runner.info(f"Generating {num_confs} conformers for molecule {molecule_smiles} / {molecule_formula} with seed {seed}")
    return generate_rdkit_conformers(molecule, num_confs=num_confs, seed=seed, profile=profile)


@node
def conformers_rdkit_from_smiles(smiles: StringData, confgen_input: ConformerGenerationInput, **kwargs) -> MoleculeList:
    node_runner = kwargs["node_runner"]
    num_confs = confgen_input.num_confs
    seed = confgen_input.seed
    profile = confgen_input.profile
    from molecular_qm_util.rdkit_scripts.smiles_to_molecule import smiles_to_molecule
    molecule = smiles_to_molecule(smiles.value)
    molecule_smiles = smiles.value
    molecule_formula = "NA"
    node_runner.info(f"Generating {num_confs} conformers for molecule {molecule_smiles} / {molecule_formula} with seed {seed}")
    return generate_rdkit_conformers(molecule, num_confs=num_confs, seed=seed, profile=profile)


def generate_rdkit_conformers(
    molecule: Molecule,
    num_confs: int = 50,
    prune_rms_thresh: float = 0.5,
    seed: int = 1,
    forcefield: str = "mmff",
    max_iters: int = 500,
    threads: int = 0,
    profile: bool = False,
) -> MoleculeList:
    """
    Generate, minimize, and rank conformers for a qm_models ``Molecule`` using RDKit.

    Returns a ``MoleculeList`` of pruned, energy-ranked conformers. Each returned
    molecule stores its force-field energy (kcal/mol) in ``properties["energy"]``
    and its energy ordering in ``properties["rank_id"]`` (set by
    :func:`molecular_qm_models.prune_conformers.prune_conformers`).
    """
    timing = {}
    t0 = time.perf_counter()

    _require_rdkit()

    smiles = molecule.smiles or molecule.make_smiles()
    print(f"Generating {num_confs} conformers for {smiles!r} using RDKit")
    logger.info(f"Using {forcefield} forcefield")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles!r}")

    mol = Chem.AddHs(mol)

    params = AllChem.ETKDGv3()
    params.randomSeed = int(seed)
    params.pruneRmsThresh = float(prune_rms_thresh) if prune_rms_thresh > 0 else -1.0
    params.useExpTorsionAnglePrefs = True
    params.useBasicKnowledge = True
    if threads > 0:
        params.numThreads = int(threads)

    # Build one initial 3D geometry separately and keep it out of mol for now.
    seed_mol = Chem.Mol(mol)
    first_id = AllChem.EmbedMolecule(seed_mol, params)
    if first_id == -1:
        # Retry with random coordinates if the first attempt fails
        print(f"RDKit embedding failed for {smiles!r}. Attempting with random coordinates...")
        retry_params = AllChem.ETKDGv3()
        retry_params.randomSeed = -1
        retry_params.useRandomCoords = True
        retry_params.pruneRmsThresh = float(prune_rms_thresh) if prune_rms_thresh > 0 else -1.0
        retry_params.useExpTorsionAnglePrefs = True
        retry_params.useBasicKnowledge = True
        if threads > 0:
            retry_params.numThreads = int(threads)

        seed_mol = Chem.Mol(mol)
        first_id = AllChem.EmbedMolecule(seed_mol, retry_params)
        if first_id == -1:
            raise RuntimeError(f"Failed to generate a 3D geometry for {smiles!r} using RDKit")
        else:
            print(f"RDKit embedding succeeded for {smiles!r}")

    # Add additional conformers to mol
    if num_confs > 1:
        AllChem.EmbedMultipleConfs(mol, numConfs=int(num_confs) - 1, params=params)

    # Add the initial conformer as well
    mol.AddConformer(seed_mol.GetConformer(int(first_id)), assignId=True)

    timing["Embedding"] = time.perf_counter() - t0
    t1 = time.perf_counter()

    all_cids = [c.GetId() for c in mol.GetConformers()]
    print(f"Generated {len(all_cids)} conformers for {smiles!r} {all_cids!r}")

    # Parallel Minimization
    if threads == 0:
        # Use half of available cores if not specified, but at least 1
        num_workers = max(1, os.cpu_count() // 2 if os.cpu_count() else 1)
    else:
        num_workers = threads

    def minimize_conformer(cid):
        # We need a local copy of mol to avoid issues if minimization is not thread-safe
        # (though RDKit forcefields usually are if they work on different conformers of the same mol,
        # but let's be safe and use a copy or just be careful)
        # Actually RDKit's MMFF/UFF minimization on DIFFERENT conformers of the SAME mol object
        # MIGHT have issues with shared state in the mol object itself if properties are being set.
        # But we are just minimizing coordinates.

        # To be absolutely safe in parallel threads with RDKit, often it's better to give each thread its own mol copy
        local_mol = Chem.Mol(mol)

        if forcefield == "mmff":
            mmff_props = AllChem.MMFFGetMoleculeProperties(local_mol, mmffVariant="MMFF94")
            ff = AllChem.MMFFGetMoleculeForceField(local_mol, mmff_props, confId=int(cid))
        else:
            ff = AllChem.UFFGetMoleculeForceField(local_mol, confId=int(cid))

        if ff is None:
            return cid, None, None

        ff.Initialize()
        ff.Minimize(maxIts=int(max_iters))
        e = float(ff.CalcEnergy())
        # Get minimized coordinates for this conformer
        conf = local_mol.GetConformer(int(cid))
        return cid, e, conf

    results = []
    print(f"Minimizing {len(all_cids)} conformers using {num_workers} threads")
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        future_to_cid = {executor.submit(minimize_conformer, cid): cid for cid in all_cids}
        for future in concurrent.futures.as_completed(future_to_cid):
            cid, e, conf = future.result()
            if e is not None:
                results.append((cid, e, conf))
                print(f"Minimized conformer {cid} for {smiles!r} - Energy: {e:.4f} kcal/mol")

    timing["Minimization"] = time.perf_counter() - t1
    t2 = time.perf_counter()

    # Sort by energy and build a qm_models Molecule for each conformer.
    results.sort(key=lambda x: x[1])

    conformers = MoleculeList()
    for cid, e, conf in results:
        conformers.append(_rdkit_conf_to_molecule(mol, conf, e, smiles))

    # Apply systematic pruning
    rms_thresh = float(prune_rms_thresh) if prune_rms_thresh > 0 else _DEFAULT_PRUNE_RMS_THRESH
    conformers = prune_conformers(conformers, rms_thresh)

    timing["Pruning"] = time.perf_counter() - t2

    if profile:
        print("\n--- RDKit CPU Profile ---")
        for k, v in timing.items():
            print(f"{k:<25}: {v:.4f}s")
        print("--------------------------\n")

    return conformers


def generate_rdkit_conformers_from_smiles(
    smiles: str,
    num_confs: int = 50,
    prune_rms_thresh: float = 0.5,
    seed: int = 1,
    forcefield: str = "mmff",
    max_iters: int = 500,
    threads: int = 0,
    profile: bool = False,
) -> MoleculeList:
    """
    SMILES wrapper around :func:`generate_rdkit_conformers`.

    Converts ``smiles`` into a qm_models ``Molecule`` and then runs the RDKit
    conformer generation on it.
    """
    from molecular_qm_util.rdkit_scripts.smiles_to_molecule import smiles_to_molecule

    molecule = smiles_to_molecule(smiles)
    molecule.smiles = smiles
    return generate_rdkit_conformers(
        molecule,
        num_confs=num_confs,
        prune_rms_thresh=prune_rms_thresh,
        seed=seed,
        forcefield=forcefield,
        max_iters=max_iters,
        threads=threads,
        profile=profile,
    )


def _write_xyz_from_rdkit(mol: "Chem.Mol", conf_id: int, out_xyz: Path) -> None:
    _require_rdkit()
    conf = mol.GetConformer(int(conf_id))

    lines = []
    lines.append(str(mol.GetNumAtoms()))
    lines.append("generated by RDKit")
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        lines.append(f"{atom.GetSymbol():<2} {pos.x: .8f} {pos.y: .8f} {pos.z: .8f}")

    out_xyz.write_text("\n".join(lines) + "\n", encoding="utf-8")
