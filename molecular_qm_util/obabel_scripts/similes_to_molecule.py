import asyncio
import os
from pathlib import Path

from simstack.models import StringData, Parameters
from simstack.core.node import node
from molecular_qm_models import Molecule, Atom
from simstack.core.context import context

try:
    from openbabel import openbabel as ob
except ImportError:
    ob = None

def _parse_xyz(xyz_text: str, *, smiles: str) -> Molecule:
    """
    Parse XYZ text into a Molecule/Atom list.

    XYZ format:
      line 1: atom count (optional-ish in the wild, but Open Babel writes it)
      line 2: comment
      line 3+: "Element  x  y  z"
    """
    xyz_text += "\n"
    lines = [ln.strip() for ln in xyz_text.splitlines() if ln.strip()]
    if len(lines) < 3:
        raise RuntimeError(f"Open Babel produced invalid/empty XYZ output for SMILES '{smiles}'")

    # Open Babel writes atom count on the first line; try to honor it but be tolerant.
    start_idx = 1
    try:
        n_atoms = int(lines[0])
        if n_atoms <= 0:
            raise ValueError
    except Exception:
        # If the first line isn't an int, treat the whole file as atom lines
        n_atoms = None
        start_idx = 0

    atom_lines = lines[start_idx:]
    if n_atoms is not None and len(atom_lines) < n_atoms:
        raise RuntimeError(
            f"Open Babel XYZ output truncated for SMILES '{smiles}': "
            f"expected {n_atoms} atoms, got {len(atom_lines)} lines"
        )

    new_molecule = Molecule()
    new_molecule.smiles = smiles

    count = n_atoms if n_atoms is not None else len(atom_lines)
    for i, ln in enumerate(atom_lines[:count], start=1):
        parts = ln.split()
        if len(parts) < 4:
            raise RuntimeError(f"Invalid XYZ atom line #{i} for SMILES '{smiles}': '{ln}'")

        element = parts[0]
        try:
            x, y, z = map(float, parts[1:4])
        except Exception as e:
            raise RuntimeError(f"Invalid XYZ coordinates on line #{i} for SMILES '{smiles}': '{ln}'") from e

        new_molecule.add_atom(Atom(element=element, x=x, y=y, z=z))

    return new_molecule

def _setup_ob_env():
    """Ensure BABEL_DATADIR is set so Open Babel can find its data files."""
    if os.environ.get("BABEL_DATADIR"):
        return

    if ob is not None:
        # Try to find data dir relative to the installed openbabel package
        try:
            import openbabel
            pkg_path = Path(openbabel.__file__).parent
            # Common path in wheels: site-packages/openbabel/bin/data
            data_dir = pkg_path / "bin" / "data"
            if data_dir.exists():
                os.environ["BABEL_DATADIR"] = str(data_dir)
                return
            
            # Another common path: site-packages/openbabel/data
            data_dir = pkg_path / "data"
            if data_dir.exists():
                os.environ["BABEL_DATADIR"] = str(data_dir)
                return
        except Exception:
            pass

@node(parameters=Parameters(resource="int-nano"), force_rerun=True)
def smiles_to_molecule_obabel(smiles: StringData, **kwargs) -> Molecule:
    """
    Converts a SMILES string into a 3D Molecule object using Open Babel's
    Python bindings.
    """
    node_runner = kwargs.get("node_runner", None)
    if node_runner:
        node_runner.info(f"Converting SMILES '{smiles}' to 3D coordinates")
    
    if isinstance(smiles, StringData):
        smiles = smiles.value
    if not smiles or not str(smiles).strip():
        raise ValueError("SMILES string cannot be empty")

    smiles = str(smiles).strip()

    if ob is None:
        raise ImportError(
            "Open Babel Python bindings not found. "
            "Please install them (e.g., pip install openbabel-wheel)."
        )

    _setup_ob_env()

    # Create OBMol and OBConversion
    mol = ob.OBMol()
    conv = ob.OBConversion()
    if not conv.SetInAndOutFormats("smi", "xyz"):
        raise RuntimeError("Open Babel: SMILES or XYZ format not supported")

    if not conv.ReadString(mol, smiles):
        raise RuntimeError(f"Open Babel failed to read SMILES: {smiles}")

    # Add Hydrogens
    mol.AddHydrogens()

    # Generate 3D coordinates
    builder = ob.OBBuilder()
    if not builder.Build(mol):
        raise RuntimeError(f"Open Babel failed to generate 3D coordinates for SMILES: {smiles}")

    # Perform forcefield minimization
    ff = ob.OBForceField.FindForceField("MMFF94")
    if ff is None:
        # Fallback to UFF if MMFF94 is missing
        ff = ob.OBForceField.FindForceField("UFF")
    
    if ff is not None:
        if ff.Setup(mol):
            ff.SteepestDescent(50)
            ff.GetCoordinates(mol)
        else:
            if node_runner:
                node_runner.warning(f"Forcefield setup failed for SMILES '{smiles}', using 3D builder coords only.")

    # Convert to XYZ text to parse via existing helper
    xyz_text = conv.WriteString(mol)
    if not xyz_text:
        raise RuntimeError(f"Open Babel produced no XYZ output for SMILES '{smiles}'")

    return _parse_xyz(xyz_text, smiles=smiles)

async def main():
    await context.initialize()
    smiles = StringData(value="C1=CC=CC=C1")
    mol = smiles_to_molecule_obabel(smiles)
    print(mol)

if __name__ == "__main__":
    asyncio.run(main())
