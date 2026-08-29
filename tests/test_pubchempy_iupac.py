from pathlib import Path
from unittest.mock import patch
import pytest
import pubchempy as pcp
from molecular_qm_models import Atom, Molecule
from molecular_qm_util.pubchempy_scripts.compute_iupac_name import compute_iupac_name


DATA_DIR = Path(__file__).parent / "data"


def test_pubchempy_get_compounds_iupac_name():
    """Test direct PubChemPy querying for IUPAC name."""
    compounds = pcp.get_compounds("c1ccccc1", namespace="smiles")
    assert compounds
    assert compounds[0].iupac_name.lower() == "benzene"


def test_compute_iupac_name_benzene_from_smi_file():
    """Test computing IUPAC name for benzene loaded from tests/data/benzene.smi."""
    smi_path = DATA_DIR / "benzene.smi"
    assert smi_path.exists()
    smiles = smi_path.read_text().strip()

    molecule = Molecule(
        smiles=smiles,
        atoms=[
            Atom(element="C", x=0.0, y=0.0, z=0.0),
            Atom(element="C", x=1.0, y=0.0, z=0.0),
            Atom(element="C", x=1.5, y=1.0, z=0.0),
            Atom(element="C", x=1.0, y=2.0, z=0.0),
            Atom(element="C", x=0.0, y=2.0, z=0.0),
            Atom(element="C", x=-0.5, y=1.0, z=0.0),
            Atom(element="H", x=0.0, y=-1.0, z=0.0),
            Atom(element="H", x=2.0, y=0.0, z=0.0),
            Atom(element="H", x=2.5, y=1.0, z=0.0),
            Atom(element="H", x=1.0, y=3.0, z=0.0),
            Atom(element="H", x=-0.5, y=3.0, z=0.0),
            Atom(element="H", x=-1.5, y=1.0, z=0.0),
        ],
    )

    iupac_name = compute_iupac_name(molecule)
    assert iupac_name.lower() == "benzene"


def test_compute_iupac_name_ethanol():
    """Test computing IUPAC name for ethanol."""
    molecule = Molecule(
        smiles="CCO",
        atoms=[
            Atom(element="C", x=0.0, y=0.0, z=0.0),
            Atom(element="C", x=1.5, y=0.0, z=0.0),
            Atom(element="O", x=2.0, y=1.2, z=0.0),
        ],
    )

    iupac_name = compute_iupac_name(molecule)
    assert iupac_name.lower() == "ethanol"


def test_compute_iupac_name_empty_molecule_raises_value_error():
    """Test that empty molecule raises ValueError."""
    empty_mol = Molecule(atoms=[])
    with pytest.raises(ValueError, match="Cannot compute IUPAC name for empty molecule"):
        compute_iupac_name(empty_mol)


def test_compute_iupac_name_fallback_formula_when_no_compounds_found():
    """Test fallback to sorted molecular formula when PubChem finds no compounds."""
    molecule = Molecule(
        smiles="UNKNOWN_SMILES_XYZ",
        atoms=[
            Atom(element="C", x=0.0, y=0.0, z=0.0),
            Atom(element="C", x=1.0, y=0.0, z=0.0),
            Atom(element="H", x=2.0, y=0.0, z=0.0),
            Atom(element="O", x=3.0, y=0.0, z=0.0),
        ],
    )

    with patch("pubchempy.get_compounds", return_value=[]):
        name = compute_iupac_name(molecule)
        assert name == "C2H1O1"


def test_compute_iupac_name_fallback_formula_on_exception():
    """Test fallback to sorted molecular formula when PubChem query raises an exception."""
    molecule = Molecule(
        smiles="c1ccccc1",
        atoms=[
            Atom(element="C", x=0.0, y=0.0, z=0.0),
            Atom(element="C", x=1.0, y=0.0, z=0.0),
            Atom(element="H", x=2.0, y=0.0, z=0.0),
        ],
    )

    with patch("pubchempy.get_compounds", side_effect=Exception("PubChem connection failed")):
        name = compute_iupac_name(molecule)
        assert name == "C2H1"


def test_molecule_make_formula():
    """Test Molecule.make_formula() calls compute_iupac_name correctly."""
    molecule = Molecule(
        smiles="c1ccccc1",
        atoms=[
            Atom(element="C", x=0.0, y=0.0, z=0.0),
            Atom(element="H", x=1.0, y=0.0, z=0.0),
        ],
    )
    formula = molecule.make_formula()
    assert formula.lower() == "benzene"
    assert molecule.formula.lower() == "benzene"
