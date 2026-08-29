import os
import shutil
import subprocess
import tempfile

import time
from pathlib import Path
from typing import Optional
from simstack.core.node import node
from molecular_qm_models import Atom, Molecule, MoleculeList
from molecular_qm_models.multi_molecule_text_parser import iter_sdf_frames
from molecular_qm_models.prune_conformers import prune_conformers
from odmantic import Model
from simstack.models import simstack_model

try:
    from openbabel import openbabel as ob
    from openbabel import pybel
except ImportError as e:  # pragma: no cover - optional dependency
    ob = None
    pybel = None
    _OPENBABEL_IMPORT_ERROR = e
else:
    _OPENBABEL_IMPORT_ERROR = None

# Default RMSD threshold (Angstrom) used when pruning OpenBabel conformers.
_DEFAULT_PRUNE_RMS_THRESH = 0.1


def _get_obabel_exe() -> Optional[str]:
    # Common locations on Windows
    paths = [
        r"C:\Program Files\Avogadro2\bin\obabel.exe",
        r"C:\Program Files (x86)\OpenBabel-3.1.1\obabel.exe",
        r"C:\Program Files\OpenBabel-3.1.1\obabel.exe",
    ]
    for p in paths:
        if Path(p).exists():
            return p
    return shutil.which("obabel")


def _require_openbabel() -> None:
    if (ob is None or pybel is None) and _get_obabel_exe() is None:
        raise RuntimeError(
            "OpenBabel is not available in this environment. Install it first "
            "(e.g. `conda install openbabel`) or ensure `obabel` is on your PATH."
        ) from _OPENBABEL_IMPORT_ERROR


def _molecule_to_pybel(molecule: Molecule):
    """
    Build a ``pybel.Molecule`` (with perceived bonds) from a qm_models Molecule.

    The molecule is serialised to an XYZ block and read back through OpenBabel so
    that bonds are perceived from the 3D geometry.
    """
    xyz_lines = [str(len(molecule.atoms)), molecule.smiles or ""]
    for atom in molecule.atoms:
        xyz_lines.append(f"{atom.element} {atom.x} {atom.y} {atom.z}")
    xyz_text = "\n".join(xyz_lines) + "\n"
    return pybel.readstring("xyz", xyz_text)


def _obmol_to_molecule(obmol, energy: float, smiles: Optional[str]) -> Molecule:
    """Convert a single OpenBabel ``OBMol`` conformer into a qm_models Molecule."""
    new_molecule = Molecule()
    new_molecule.smiles = smiles
    for atom in ob.OBMolAtomIter(obmol):
        new_molecule.add_atom(
            Atom(
                element=ob.GetSymbol(atom.GetAtomicNum()),
                x=atom.GetX(),
                y=atom.GetY(),
                z=atom.GetZ(),
            )
        )
    new_molecule.properties["energy"] = energy
    return new_molecule


def _energy_from_molblock(molblock: str) -> float:
    """Extract the ``Energy`` SDF tag from a molblock, defaulting to ``0.0``."""
    lines = molblock.splitlines()
    for idx, line in enumerate(lines):
        if line.strip().lower().startswith("> <energy"):
            if idx + 1 < len(lines):
                try:
                    return float(lines[idx + 1].strip())
                except ValueError:
                    return 0.0
    return 0.0


@simstack_model
class ConformerGenerationInput(Model):
    """Input for :func:`generate_openbabel_conformers`."""
    num_confs: int = 50
    seed: int = 1
    profile: bool = False

@node
def conformers_openbabel(molecule: Molecule, confgen_input: ConformerGenerationInput,**kwargs) -> MoleculeList:
    node_runner = kwargs["node_runner"]
    num_confs = confgen_input.num_confs
    seed = confgen_input.seed
    profile = confgen_input.profile
    molecule_smiles = molecule.smiles if molecule.smiles is not None else "NA"
    molecule_formula = molecule.formula if molecule.formula is not None else "NA"
    node_runner.info(f"Generating {num_confs} conformers for molecule {molecule_smiles} / {molecule_formula} with seed {seed}")
    return generate_openbabel_conformers(molecule,num_confs=num_confs,seed=seed,profile=profile)

@node
def conformers_openbabel_from_smiles(smiles: str, confgen_input: ConformerGenerationInput,**kwargs) -> MoleculeList:
    node_runner = kwargs["node_runner"]
    num_confs = confgen_input.num_confs
    seed = confgen_input.seed
    profile = confgen_input.profile
    from molecular_qm_util.rdkit_scripts.smiles_to_molecule import smiles_to_molecule
    molecule = smiles_to_molecule(smiles)
    molecule_smiles = smiles
    molecule_formula = "NA"


    node_runner.info(f"Generating {num_confs} conformers for molecule {molecule_smiles} / {molecule_formula} with seed {seed}")
    return generate_openbabel_conformers(molecule,num_confs=num_confs,seed=seed,profile=profile)



def generate_openbabel_conformers(
    molecule: Molecule,
    num_confs: int = 50,
    seed: int = 1,
    profile: bool = False,
) -> MoleculeList:
    """
    Generate 3D conformers for a qm_models ``Molecule`` using OpenBabel.

    Returns a ``MoleculeList`` of pruned, energy-ranked conformers. Each returned
    molecule stores its force-field energy (kcal/mol) in ``properties["energy"]``
    and its energy ordering in ``properties["rank_id"]`` (set by
    :func:`molecular_qm_models.prune_conformers.prune_conformers`).
    """
    _require_openbabel()

    timing: dict = {}

    if ob is not None and pybel is not None:
        ranked = _generate_openbabel_conformers_python(molecule, num_confs, seed, timing)
    else:
        ranked = _generate_openbabel_conformers_cli(molecule, num_confs, seed, timing)

    if profile:
        print("\n--- OpenBabel CPU Profile ---")
        for k, v in timing.items():
            print(f"{k:<25}: {v:.4f}s")
        print("-----------------------------\n")

    return ranked


def generate_openbabel_conformers_from_smiles(
    smiles: str,
    num_confs: int = 50,
    seed: int = 1,
    profile: bool = False,
) -> MoleculeList:
    """
    SMILES wrapper around :func:`generate_openbabel_conformers`.

    Converts ``smiles`` into a qm_models ``Molecule`` and then runs the
    OpenBabel conformer generation on it.
    """
    from molecular_qm_util.rdkit_scripts.smiles_to_molecule import smiles_to_molecule

    molecule = smiles_to_molecule(smiles)
    molecule.smiles = smiles
    return generate_openbabel_conformers(
        molecule, num_confs=num_confs, seed=seed, profile=profile
    )


def _generate_openbabel_conformers_python(
    molecule: Molecule,
    num_confs: int = 50,
    seed: int = 1,
    timing: dict = None,
) -> MoleculeList:
    t0 = time.perf_counter()
    smiles = molecule.smiles
    mol = _molecule_to_pybel(molecule)
    mol.addh()

    # Find forcefield for initial minimization and final ranking
    ff = ob.OBForceField.FindForceField("mmff94")
    if not ff:
        ff = ob.OBForceField.FindForceField("uff")

    # Set random seed if provided
    if seed:
        ob.obErrorLog.StopLogging()
        # OpenBabel's OBConformerSearch doesn't have a direct SetSeed,
        # but we can try to seed the global RNG if needed.
        # However, OBConformerSearch is somewhat deterministic by default or uses its own logic.

    # 1. Minimize the initial conformer to use as a seed
    if ff:
        success = ff.Setup(mol.OBMol)
        if success:
            ff.ConjugateGradients(500)
            ff.GetCoordinates(mol.OBMol)
            print(f"Initial seed minimization for {smiles!r} - Energy: {ff.Energy():.4f} kcal/mol")

    if timing is not None:
        timing["Initial setup/Min"] = time.perf_counter() - t0
    t1 = time.perf_counter()

    # Use OBConformerSearch for multiple conformers
    if num_confs > 1:
        # Stop logging to avoid cluttering stderr with OBConformerSearch warnings
        ob.obErrorLog.StopLogging()
        cs = ob.OBConformerSearch()
        # Setup(OBMol mol, int numConformers, int numChildren, int mutability, int convergence)
        cs.Setup(mol.OBMol, num_confs, 50, 5, 25)
        cs.Search()
        cs.GetConformers(mol.OBMol)
        # Restart logging
        ob.obErrorLog.StartLogging()

        # Check if OBConformerSearch failed to generate any conformers despite search
        if mol.OBMol.NumConformers() == 0:
            print(f"OpenBabel's OBConformerSearch failed for {smiles!r}. Using initial 3D structure.")

    if timing is not None:
        timing["Conformer Search"] = time.perf_counter() - t1
    t2 = time.perf_counter()

    conformers = MoleculeList()
    num_found = mol.OBMol.NumConformers()
    if num_found == 0:
        # If no conformers found by Search, use the initial 3D one
        num_found = 1

    # Setup for each conformer as minimized separately
    for i in range(num_found):
        mol.OBMol.SetConformer(i)
        energy = 0.0
        if ff:
            success = ff.Setup(mol.OBMol)
            if success:
                # Use fewer steps or check for NaN
                ff.ConjugateGradients(500)
                ff.GetCoordinates(mol.OBMol)
                energy = ff.Energy()
                print(f"Minimized OpenBabel conformer {i} for {smiles!r} - Energy: {energy:.4f} kcal/mol")
                if float('inf') == energy or energy != energy:  # Check for Inf or NaN
                    energy = 1e10  # Use a large value for failed minimizations
            else:
                # If Setup fails, maybe use basic coords and a very high energy
                energy = 1e10

        # Create a copy of the molecule for each conformer to store the coordinates.
        # Use ob.OBMol(mol.OBMol) to ensure a deep copy of the structure.
        new_ob_mol = ob.OBMol(mol.OBMol)
        new_ob_mol.SetConformer(i)
        conformers.append(_obmol_to_molecule(new_ob_mol, energy, smiles))

    if timing is not None:
        timing["Final Minimization"] = time.perf_counter() - t2
    t3 = time.perf_counter()

    print(f"Generated {len(conformers)} conformers for {smiles!r} using OpenBabel Python bindings")
    conformers = prune_conformers(conformers, _DEFAULT_PRUNE_RMS_THRESH)

    if timing is not None:
        timing["Pruning"] = time.perf_counter() - t3

    return conformers


def _generate_openbabel_conformers_cli(
    molecule: Molecule,
    num_confs: int = 50,
    seed: int = 1,
    timing: dict = None,
) -> MoleculeList:
    t0 = time.perf_counter()
    ob_exe = _get_obabel_exe()
    if not ob_exe:
        raise RuntimeError("OpenBabel executable not found.")

    smiles = molecule.smiles

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_xyz = os.path.join(tmpdir, "input.xyz")
        molecule.to_file(tmp_xyz)

        tmp_sdf = os.path.join(tmpdir, "conformers.sdf")

        # Use obabel CLI to generate conformers.
        # --confsearch performs a conformer search
        # --nconf specifies the number of conformers
        cmd = [
            ob_exe,
            tmp_xyz,
            "-O", tmp_sdf,
            "--confsearch",
            "--nconf", str(num_confs),
            "--weighted",  # Optional: weight by energy
            "--ff", "mmff94",
        ]

        print(f"Running OpenBabel CLI: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"OpenBabel CLI failed: {result.stderr}")
            # Try without confsearch if it fails, just minimize the input geometry
            cmd = [ob_exe, tmp_xyz, "-O", tmp_sdf, "--ff", "mmff94"]
            subprocess.run(cmd, capture_output=True, text=True)

        if timing is not None:
            timing["CLI Execution"] = time.perf_counter() - t0
        t1 = time.perf_counter()

        conformers = MoleculeList()
        if os.path.exists(tmp_sdf):
            with open(tmp_sdf, "r") as f:
                sdf_text = f.read()
            for molblock in iter_sdf_frames(sdf_text):
                conf = Molecule.from_sdf(molblock)
                conf.smiles = smiles
                conf.properties["energy"] = _energy_from_molblock(molblock)
                conformers.append(conf)

        print(f"Generated {len(conformers)} conformers for {smiles!r} using OpenBabel CLI")
        conformers = prune_conformers(conformers, _DEFAULT_PRUNE_RMS_THRESH)

        if timing is not None:
            timing["Parsing/Pruning"] = time.perf_counter() - t1

        return conformers
