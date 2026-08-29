# ``compute_iupac_name`` only depends on pubchempy (and molecular_qm_models),
# so it is imported unconditionally. The obabel- and rdkit-based helpers pull in
# heavy optional dependencies (openbabel, rdkit, simstack); guard them so that a
# missing optional dependency does not break importing the rest of the package.
from .pubchempy_scripts.compute_iupac_name import compute_iupac_name

__all__ = ["compute_iupac_name"]

try:
    from .obabel_scripts.similes_to_molecule import smiles_to_molecule
except ImportError:
    pass
else:
    __all__.append("smiles_to_molecule")

try:
    from .rdkit_scripts.molecule_to_smiles import (
        compute_smiles,
        molecule_to_canonical_smiles,
        molecule_to_smiles,
    )
except ImportError:
    pass
else:
    __all__ += [
        "compute_smiles",
        "molecule_to_smiles",
        "molecule_to_canonical_smiles",
    ]
