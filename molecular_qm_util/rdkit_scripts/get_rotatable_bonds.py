from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

def get_rotatable_bonds(mol: Chem.Mol) -> List[Tuple[int, int, int, int]]:
    """
    Find rotatable bonds (dihedral angles) in a molecule using RDKit.
    Returns a list of 4-atom indices representing the dihedrals.
    """
    rotatable_bond_smarts = "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"
    rotatable_bond_query = Chem.MolFromSmarts(rotatable_bond_smarts)
    matches = mol.GetSubstructMatches(rotatable_bond_query)

    dihedrals = []
    for atom1_idx, atom2_idx in matches:
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)
        n1 = [neighbor.GetIdx() for neighbor in atom1.GetNeighbors() if neighbor.GetIdx() != atom2_idx]
        n2 = [neighbor.GetIdx() for neighbor in atom2.GetNeighbors() if neighbor.GetIdx() != atom1_idx]
        if n1 and n2:
            dihedrals.append((n1[0], atom1_idx, atom2_idx, n2[0]))
    return dihedrals

