from rdkit import Chem
from rdkit.Chem import AllChem
from molecular_qm_models import Molecule


def smiles_to_molecule(smiles_string: str) -> Molecule:
    mol = Chem.MolFromSmiles(smiles_string)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    return Molecule.from_sdf(Chem.MolToMolBlock(mol))
