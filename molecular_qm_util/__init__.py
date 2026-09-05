from .pubchempy_scripts.compute_iupac_name import compute_iupac_name
from .rdkit_scripts.smiles_to_molecule import smiles_to_molecule
from .rdkit_scripts.molecule_to_smiles import (
    compute_smiles,
    molecule_to_canonical_smiles,
    molecule_to_smiles,
)
from .rdkit_scripts.rdkit_optimize import (
    rdkit_optimize,
    RDKitForceField,
    molecule_to_rdkit,
    rdkit_to_molecule,
)
from .rdkit_scripts.get_rotatable_bonds import (
    get_rotatable_bonds,
    get_rotatable_bonds_base,
    get_rotable_bonds,
)

__all__ = [
    "compute_iupac_name",
    "smiles_to_molecule",
    "compute_smiles",
    "molecule_to_smiles",
    "molecule_to_canonical_smiles",
    "rdkit_optimize",
    "RDKitForceField",
    "molecule_to_rdkit",
    "rdkit_to_molecule",
    "get_rotatable_bonds",
    "get_rotatable_bonds_base",
    "get_rotable_bonds",
]

try:
    from .obabel_scripts.similes_to_molecule import smiles_to_molecule_obabel
except ImportError:
    pass
else:
    __all__.append("smiles_to_molecule_obabel")
