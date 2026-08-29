from .pubchempy_scripts.compute_iupac_name import compute_iupac_name
from .rdkit_scripts.smiles_to_molecule import smiles_to_molecule
from .rdkit_scripts.molecule_to_smiles import (
    compute_smiles,
    molecule_to_canonical_smiles,
    molecule_to_smiles,
)

__all__ = [
    "compute_iupac_name",
    "smiles_to_molecule",
    "compute_smiles",
    "molecule_to_smiles",
    "molecule_to_canonical_smiles",
]

try:
    from .obabel_scripts.similes_to_molecule import smiles_to_molecule_obabel
except ImportError:
    pass
else:
    __all__.append("smiles_to_molecule_obabel")
