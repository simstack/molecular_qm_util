import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from molecular_qm_models import Atom, Molecule
from molecular_qm_util.rdkit_scripts import (
    compute_smiles,
    molecule_to_canonical_smiles,
    molecule_to_smiles,
)
from molecular_qm_util.rdkit_scripts.compute_smiles import (
    compute_smiles as compute_smiles_from_module,
)
from molecular_qm_util.rdkit_scripts.molecule_to_smiles import (
    molecule_to_smiles as molecule_to_smiles_from_module,
)
from molecular_qm_util import compute_smiles as top_compute_smiles


def _make_3d_molecule(smiles: str, charge: int = 0) -> Molecule:
    rd_mol = Chem.MolFromSmiles(smiles)
    assert rd_mol is not None
    rd_mol = Chem.AddHs(rd_mol)
    AllChem.EmbedMolecule(rd_mol, AllChem.ETKDG())
    conf = rd_mol.GetConformer()
    atoms = []
    for atom in rd_mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append(
            Atom(
                element=atom.GetSymbol(),
                x=float(pos.x),
                y=float(pos.y),
                z=float(pos.z),
            )
        )
    return Molecule(atoms=atoms, properties={"charge": charge})


@pytest.mark.parametrize(
    "smiles",
    [
        "CCO",
        "c1ccccc1",
        "CC(=O)O",
        "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        "CC#N",
        "C#C",
        "C=C",
        "O",
        "N",
        "C",
    ],
)
def test_molecule_to_canonical_smiles(smiles: str):
    molecule = _make_3d_molecule(smiles)
    canon_smiles = molecule_to_smiles(molecule)
    expected = Chem.CanonSmiles(smiles)
    assert canon_smiles == expected


def test_alias_functions():
    molecule = _make_3d_molecule("CCO")
    assert compute_smiles(molecule) == "CCO"
    assert molecule_to_canonical_smiles(molecule) == "CCO"
    assert compute_smiles_from_module(molecule) == "CCO"
    assert molecule_to_smiles_from_module(molecule) == "CCO"
    assert top_compute_smiles(molecule) == "CCO"


def test_molecule_make_smiles_integration():
    molecule = _make_3d_molecule("c1ccccc1")
    smiles = molecule.make_smiles()
    assert smiles == "c1ccccc1"
    assert molecule.smiles == "c1ccccc1"


def test_chiral_molecules():
    l_alanine = _make_3d_molecule("N[C@@H](C)C(=O)O")
    d_alanine = _make_3d_molecule("N[C@H](C)C(=O)O")
    assert molecule_to_smiles(l_alanine) == Chem.CanonSmiles("N[C@@H](C)C(=O)O")
    assert molecule_to_smiles(d_alanine) == Chem.CanonSmiles("N[C@H](C)C(=O)O")
    assert molecule_to_smiles(l_alanine) != molecule_to_smiles(d_alanine)


def test_charged_molecules():
    nh4 = _make_3d_molecule("[NH4+]", charge=1)
    assert molecule_to_smiles(nh4) == "[NH4+]"

    acetate = _make_3d_molecule("CC(=O)[O-]", charge=-1)
    assert molecule_to_smiles(acetate) == "CC(=O)[O-]"


def test_empty_molecule_raises():
    empty_mol = Molecule()
    with pytest.raises(ValueError, match="Cannot compute SMILES for empty molecule"):
        molecule_to_smiles(empty_mol)

    none_mol = None
    with pytest.raises(ValueError, match="Cannot compute SMILES for empty molecule"):
        molecule_to_smiles(none_mol)


def test_invalid_atom_raises():
    invalid_mol = Molecule(atoms=[Atom(element="", x=0.0, y=0.0, z=0.0)])
    with pytest.raises(ValueError, match="Invalid atom"):
        molecule_to_smiles(invalid_mol)
