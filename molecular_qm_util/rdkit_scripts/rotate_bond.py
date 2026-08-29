import sys
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms


def get_fragment_atoms(mol, bond_idx, start_atom_idx):
    """Finds all atoms on one side of a bond."""

    atoms_to_move = []
    stack = [start_atom_idx]
    visited = {mol.GetBondWithIdx(bond_idx).GetBeginAtomIdx(), 
               mol.GetBondWithIdx(bond_idx).GetEndAtomIdx()}
    
    while stack:
        curr = stack.pop()
        if curr not in atoms_to_move:
            atoms_to_move.append(curr)
            for neighbor in mol.GetAtomWithIdx(curr).GetNeighbors():
                idx = neighbor.GetIdx()
                if idx not in visited:
                    stack.append(idx)
                    visited.add(idx)
    return atoms_to_move

def rotate_bond_rigidly(smiles, atom1_idx, atom2_idx, steps=12):
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    
    bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
    if not bond:
        raise ValueError(f"No bond between {atom1_idx} and {atom2_idx}")
    
    # 1. Identify which atoms are in the 'moving' fragment (side of atom2)
    moving_atoms = get_fragment_atoms(mol, bond.GetIdx(), atom2_idx)
    
    # 2. Define the rotation axis (the bond vector)
    conf = mol.GetConformer()
    p1 = np.array(conf.GetAtomPosition(atom1_idx))
    p2 = np.array(conf.GetAtomPosition(atom2_idx))
    axis = p2 - p1
    
    conformers = []
    base_coords = np.array(conf.GetPositions())
    
    for i in range(steps):
        angle_rad = (2 * np.pi / steps) * i
        
        # Create a new conformer based on the original
        new_conf = Chem.Conformer(conf)
        
        # Apply rotation matrix to the fragment
        # We translate p2 to origin, rotate around axis, then translate back
        translation = p2
        rotation = rdMolTransforms.GetRotationMatrixAroundStartAndEndPoints(p1, p2, angle_rad * 180 / np.pi)
        
        for idx in moving_atoms:
            pos = base_coords[idx]
            # Transform point: (R * (p - p1)) + p1
            new_pos = rotation.dot(pos - p1) + p1
            new_conf.SetAtomPosition(idx, new_pos)
            
        conformers.append(Chem.MolToMolBlock(mol, confId=mol.AddConformer(new_conf, assignId=True)))
        
    return conformers

if __name__ == "__main__":
    smiles_in = sys.argv[1]
    a1, a2 = int(sys.argv[2]), int(sys.argv[3])
    results = rotate_bond_rigidly(smiles_in, a1, a2)
    print("$$$$".join(results))
