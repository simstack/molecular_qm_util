import sys
from rdkit import Chem
from rdkit.Chem import AllChem

def get_3d(smiles):

    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    return Chem.MolToMolBlock(mol)

if __name__ == "__main__":
    print(get_3d(sys.argv[1]))
