import importlib

__all__ = [
    "cdx_to_molecule",
    "smiles_to_molecule_obabel",
    "ConformerGenerationInput",
    "conformers_openbabel",
    "conformers_openbabel_from_smiles",
    "generate_openbabel_conformers",
    "generate_openbabel_conformers_from_smiles",
]

_EXPORTS = {
    "cdx_to_molecule": ("cdx_to_molecule", "cdx_to_molecule"),
    "smiles_to_molecule_obabel": ("similes_to_molecule", "smiles_to_molecule_obabel"),
    "ConformerGenerationInput": ("openbabel_conformers", "ConformerGenerationInput"),
    "conformers_openbabel": ("openbabel_conformers", "conformers_openbabel"),
    "conformers_openbabel_from_smiles": (
        "openbabel_conformers",
        "conformers_openbabel_from_smiles",
    ),
    "generate_openbabel_conformers": (
        "openbabel_conformers",
        "generate_openbabel_conformers",
    ),
    "generate_openbabel_conformers_from_smiles": (
        "openbabel_conformers",
        "generate_openbabel_conformers_from_smiles",
    ),
}


def __getattr__(name: str):
    if name not in _EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    module_name, attr_name = _EXPORTS[name]
    module = importlib.import_module(f".{module_name}", __name__)
    value = getattr(module, attr_name)
    globals()[name] = value
    return value
