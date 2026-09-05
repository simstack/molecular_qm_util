from __future__ import annotations
from typing import List, Optional, Tuple, Union

try:
    from rdkit import Chem
except ImportError:
    Chem = None

try:
    from molecular_qm_models import (
        Molecule,
        InternalCoordinatesList,
        InternalDihedralCoordinate,
    )
    from molecular_qm_models.internal_coordinates import InternalCoordinateBondType
except ImportError:
    Molecule = None
    InternalCoordinatesList = None
    InternalDihedralCoordinate = None
    InternalCoordinateBondType = None

from .rdkit_optimize import molecule_to_rdkit, rdkit_to_molecule


def get_rotatable_bonds_base(
    mol: Chem.Mol,
    match_double_bonds: bool = False,
    include_bond_types: bool = False,
) -> Union[List[Tuple[int, int, int, int]], Tuple[List[Tuple[int, int, int, int]], List[str]]]:
    """
    Find rotatable bonds (dihedral angles) in a molecule using RDKit.

    :param mol: RDKit molecule
    :param match_double_bonds: Whether to match double bonds as well as single bonds
    :param include_bond_types: If True, returns a tuple (dihedrals, dihedral_types)
    :return: List of 4-atom indices representing dihedrals, or (dihedrals, dihedral_types)
    """
    if Chem is None:
        raise RuntimeError("RDKit is required to find rotatable bonds.")

    use_any_bond = match_double_bonds or include_bond_types
    rotatable_bond_smarts = "[!$(*#*)&!D1]~&!@[!$(*#*)&!D1]" if use_any_bond else "[!$(*#*)&!D1]-&!@[!$(*#*)&!D1]"
    rotatable_bond_query = Chem.MolFromSmarts(rotatable_bond_smarts)
    matches = mol.GetSubstructMatches(rotatable_bond_query)

    dihedrals = []
    dihedral_types = []

    for atom1_idx, atom2_idx in matches:
        atom1 = mol.GetAtomWithIdx(atom1_idx)
        atom2 = mol.GetAtomWithIdx(atom2_idx)

        bond = mol.GetBondBetweenAtoms(atom1_idx, atom2_idx)
        if bond is not None:
            bond_order = bond.GetBondTypeAsDouble()
            dihedral_type = 'DB' if bond_order >= 1.5 else 'SB'
        else:
            dihedral_type = 'SB'

        dihedral_types.append(dihedral_type)

        n1 = [neighbor.GetIdx() for neighbor in atom1.GetNeighbors() if neighbor.GetIdx() != atom2_idx]
        n2 = [neighbor.GetIdx() for neighbor in atom2.GetNeighbors() if neighbor.GetIdx() != atom1_idx]
        if n1 and n2:
            dihedrals.append((n1[0], atom1_idx, atom2_idx, n2[0]))

    if include_bond_types:
        return dihedrals, dihedral_types
    return dihedrals


def get_rotatable_bonds(
    molecule: Molecule,
    match_double_bonds: bool = False,
    min_value: float = -180.0,
    max_value: float = 180.0,
) -> InternalCoordinatesList:
    """
    Find rotatable bonds in a Molecule and return them as an InternalCoordinatesList.

    :param molecule: Input Molecule object from molecular_qm_models
    :param match_double_bonds: Whether to match double bonds as well as single bonds
    :param min_value: Minimum dihedral value in degrees (default -180.0)
    :param max_value: Maximum dihedral value in degrees (default 180.0)
    :return: InternalCoordinatesList containing InternalDihedralCoordinate elements
    """
    if InternalCoordinatesList is None or InternalDihedralCoordinate is None:
        raise RuntimeError("molecular_qm_models is required to create InternalCoordinatesList.")

    if isinstance(molecule, Molecule):
        rd_mol = molecule_to_rdkit(molecule)
        mol_obj = molecule
    elif Chem is not None and isinstance(molecule, Chem.Mol):
        rd_mol = molecule
        mol_obj = rdkit_to_molecule(molecule)
    else:
        raise TypeError(f"Expected a Molecule object, got {type(molecule)}")

    dihedrals, dihedral_types = get_rotatable_bonds_base(
        rd_mol,
        match_double_bonds=match_double_bonds,
        include_bond_types=True,
    )

    coordinates = []
    for (a1, a2, a3, a4), dtype in zip(dihedrals, dihedral_types):
        dc = InternalDihedralCoordinate.initialize(
            a1, a2, a3, a4,
            min_value=min_value,
            max_value=max_value,
            molecule=mol_obj,
        )
        if InternalCoordinateBondType is not None:
            dc.bond_type = (
                InternalCoordinateBondType.DOUBLE
                if dtype == 'DB'
                else InternalCoordinateBondType.SINGLE
            )
        try:
            dc.compute(mol_obj)
        except Exception:
            pass
        coordinates.append(dc)

    return InternalCoordinatesList(elements=coordinates)


# Alias for spelling variations
get_rotable_bonds = get_rotatable_bonds
