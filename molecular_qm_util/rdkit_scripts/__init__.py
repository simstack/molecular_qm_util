from .generate_3d_structure_from_smiles import get_3d
from .get_rotatable_bonds import get_rotatable_bonds, get_rotatable_bonds_base, get_rotable_bonds
from .molecule_to_smiles import (
    compute_smiles,
    molecule_to_canonical_smiles,
    molecule_to_smiles,
)
from .rotate_bond import get_fragment_atoms, rotate_bond_rigidly
from .smiles_to_molecule import smiles_to_molecule
from .rdkit_optimize import (
    rdkit_optimize,
    RDKitForceField,
    RDKitForceFieldParams,
    RDKitOptimizationParams,
    molecule_to_rdkit,
    rdkit_to_molecule,
    evaluate_molecules_rdkit,
    score_molecules_rdkit,
    optimize_molecules_rdkit,
)

__all__ = [
    "get_3d",
    "get_rotatable_bonds",
    "get_rotatable_bonds_base",
    "get_rotable_bonds",
    "get_fragment_atoms",
    "rotate_bond_rigidly",
    "smiles_to_molecule",
    "molecule_to_smiles",
    "compute_smiles",
    "molecule_to_canonical_smiles",
    "ConformerGenerationInput",
    "conformers_rdkit",
    "conformers_rdkit_from_smiles",
    "generate_rdkit_conformers",
    "generate_rdkit_conformers_from_smiles",
    "rdkit_optimize",
    "RDKitForceField",
    "RDKitForceFieldParams",
    "RDKitOptimizationParams",
    "molecule_to_rdkit",
    "rdkit_to_molecule",
    "evaluate_molecules_rdkit",
    "score_molecules_rdkit",
    "optimize_molecules_rdkit",
]

_CONFORMER_EXPORTS = {
    "ConformerGenerationInput",
    "conformers_rdkit",
    "conformers_rdkit_from_smiles",
    "generate_rdkit_conformers",
    "generate_rdkit_conformers_from_smiles",
}


def __getattr__(name: str):
    if name not in _CONFORMER_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from .rdkit_conformers import (
        ConformerGenerationInput,
        conformers_rdkit,
        conformers_rdkit_from_smiles,
        generate_rdkit_conformers,
        generate_rdkit_conformers_from_smiles,
    )

    value = {
        "ConformerGenerationInput": ConformerGenerationInput,
        "conformers_rdkit": conformers_rdkit,
        "conformers_rdkit_from_smiles": conformers_rdkit_from_smiles,
        "generate_rdkit_conformers": generate_rdkit_conformers,
        "generate_rdkit_conformers_from_smiles": generate_rdkit_conformers_from_smiles,
    }[name]
    globals()[name] = value
    return value
