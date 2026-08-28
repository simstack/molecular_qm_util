from molecular_qm_models import Atom, Molecule


def rdkit_mol_to_molecule(mol, conf_id: int) -> Molecule:
    conf = mol.GetConformer(conf_id)
    atoms = []
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        atoms.append(Atom(
            element=atom.GetSymbol(),
            x=pos.x,
            y=pos.y,
            z=pos.z
        ))
    return Molecule(atoms=atoms)


def pybel_mol_to_molecule(mol) -> Molecule:
    try:
        from openbabel import pybel
    except ImportError:
        pybel = None

    if pybel is None:
        raise RuntimeError(
            "OpenBabel Python bindings (pybel) are required for pybel molecule conversion."
        )
    atoms = []
    for atom in mol.atoms:
        pos = atom.coords
        atoms.append(Atom(
            element=pybel.element_table[atom.atomicnum],
            x=pos[0],
            y=pos[1],
            z=pos[2]
        ))
    return Molecule(atoms=atoms)


def simstack_molecule_to_rdkit(mol: Molecule):
    """Convert a SimStack Molecule to an RDKit Mol via OpenBabel/pybel.

    This builds an XYZ string from the Molecule, reads it with pybel, and then
    converts to an RDKit molecule for use as an initial structure in GA.
    """

    try:
        from openbabel import pybel as ob_pybel
    except ImportError:
        ob_pybel = None

    try:
        from rdkit import Chem
    except ImportError:
        Chem = None

    if ob_pybel is None:
        raise RuntimeError(
            "OpenBabel Python bindings (pybel) are required to use an initial Molecule in GA."
        )
    if Chem is None:
        raise RuntimeError("RDKit is required to use an initial Molecule in GA.")

    n_atoms = len(mol.atoms)
    lines = [str(n_atoms), "SimStack Molecule"]
    for atom in mol.atoms:
        lines.append(f"{atom.element} {atom.x:.8f} {atom.y:.8f} {atom.z:.8f}")
    xyz_str = "\n".join(lines)

    ob_mol = ob_pybel.readstring("xyz", xyz_str)
    ob_mol.addh()
    ob_mol.make3D()
    ob_sdf = ob_mol.write("sdf")
    rd_mol = Chem.MolFromMolBlock(ob_sdf, removeHs=False)
    if rd_mol is None:
        raise RuntimeError("Failed to convert SimStack Molecule to RDKit Mol via OpenBabel.")
    return rd_mol