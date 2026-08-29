import pytest

from molecular_qm_models import Molecule, MoleculeList

from molecular_qm_util.obabel_scripts import openbabel_conformers
from molecular_qm_util.obabel_scripts.openbabel_conformers import (
    generate_openbabel_conformers,
    generate_openbabel_conformers_from_smiles,
    _energy_from_molblock,
    _get_obabel_exe,
    _molecule_to_pybel,
    _require_openbabel,
)
from molecular_qm_util.rdkit_scripts.smiles_to_molecule import smiles_to_molecule

# Ethanol (CCO) has 2 C + 6 H + 1 O = 9 atoms once hydrogens are added.
_ETHANOL_SMILES = "CCO"
_ETHANOL_NUM_ATOMS = 9

# OpenBabel is an optional dependency; skip the tests that require it when the
# python bindings (or the CLI executable) are unavailable.
_OPENBABEL_AVAILABLE = (
    openbabel_conformers.ob is not None and openbabel_conformers.pybel is not None
) or _get_obabel_exe() is not None

_requires_openbabel = pytest.mark.skipif(
    not _OPENBABEL_AVAILABLE, reason="OpenBabel is not available in this environment"
)


def _assert_ranked_molecule_list(conformers, expected_atoms):
    """Common assertions for the ``MoleculeList`` returned by the generators."""
    assert isinstance(conformers, MoleculeList)
    assert len(conformers) >= 1

    energies = []
    for molecule in conformers:
        assert isinstance(molecule, Molecule)
        assert len(molecule.atoms) == expected_atoms
        assert "energy" in molecule.properties
        assert isinstance(molecule.properties["energy"], float)
        assert "rank_id" in molecule.properties
        energies.append(molecule.properties["energy"])

    rank_ids = [m.properties["rank_id"] for m in conformers]
    assert rank_ids == list(range(len(conformers)))
    assert energies == sorted(energies)


def test_energy_from_molblock_reads_energy_tag():
    molblock = "\n".join(
        [
            "some header",
            "M  END",
            "> <Energy>",
            "-42.5",
            "",
            "$$$$",
        ]
    )
    assert _energy_from_molblock(molblock) == -42.5


def test_energy_from_molblock_defaults_to_zero_when_missing():
    molblock = "\n".join(["some header", "M  END", "$$$$"])
    assert _energy_from_molblock(molblock) == 0.0


def test_energy_from_molblock_non_numeric_defaults_to_zero():
    molblock = "\n".join(["> <Energy>", "not_a_number", "$$$$"])
    assert _energy_from_molblock(molblock) == 0.0


def test_get_obabel_exe_returns_str_or_none():
    result = _get_obabel_exe()
    assert result is None or isinstance(result, str)


@_requires_openbabel
def test_require_openbabel_available():
    # Must not raise when OpenBabel is available.
    _require_openbabel()


@_requires_openbabel
def test_molecule_to_pybel_builds_molecule():
    if openbabel_conformers.pybel is None:
        pytest.skip("OpenBabel python bindings not available")

    molecule = smiles_to_molecule(_ETHANOL_SMILES)
    molecule.smiles = _ETHANOL_SMILES

    pybel_mol = _molecule_to_pybel(molecule)
    # Hydrogens may be implicit before addh(); heavy atoms must be present.
    assert pybel_mol.OBMol.NumAtoms() >= 3


@_requires_openbabel
def test_generate_openbabel_conformers_returns_ranked_molecule_list():
    molecule = smiles_to_molecule(_ETHANOL_SMILES)
    molecule.smiles = _ETHANOL_SMILES

    conformers = generate_openbabel_conformers(molecule, num_confs=3, seed=1)

    _assert_ranked_molecule_list(conformers, _ETHANOL_NUM_ATOMS)


@_requires_openbabel
def test_generate_openbabel_conformers_from_smiles():
    conformers = generate_openbabel_conformers_from_smiles(
        _ETHANOL_SMILES, num_confs=3, seed=1
    )

    _assert_ranked_molecule_list(conformers, _ETHANOL_NUM_ATOMS)
    for molecule in conformers:
        assert molecule.smiles == _ETHANOL_SMILES
