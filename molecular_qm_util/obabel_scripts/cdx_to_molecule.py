import asyncio
import os
from pathlib import Path

from simstack.models.files import FileStack
from simstack.models import Parameters
from simstack.core.node import node
from molecular_qm_models import Molecule, Atom
from simstack.core.context import context

try:
    from openbabel import openbabel as ob
except ImportError:
    ob = None

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

#@node(parameters=Parameters(resource="int-nano"), force_rerun=True)
def cdx_to_molecule(cdx_file: FileStack, **kwargs) -> Molecule:
    """
    Convert a ChemDraw CDX file to a Molecule object using OpenBabel.

    :param cdx_file: The CDX file (as FileStack)
    :return: A Molecule object
    """
    node_runner = kwargs.get("node_runner", None)
    
    file_name = Path(cdx_file.name)
    if file_name.exists():
        file_name.unlink()

    file_path = cdx_file.get()
    if node_runner:
        node_runner.info(f"Converting CDX file '{file_path}' to Molecule")

    if ob is None:
        raise ImportError(
            "Open Babel Python bindings not found. "
            "Please install them (e.g., pip install openbabel-wheel)."
        )

    _setup_ob_env()

    try:
        conv = ob.OBConversion()
        if not conv.SetInFormat("cdx"):
            raise RuntimeError("Open Babel: CDX format not supported")
            
        mol = ob.OBMol()
        if not conv.ReadFile(mol, str(file_path)):
            raise RuntimeError(f"Open Babel failed to read CDX file: {file_path}")

        new_molecule = Molecule()
        
        for atom in ob.OBMolAtomIter(mol):
            new_molecule.add_atom(Atom(
                element=atom.GetAtomicSymbol(),
                x=atom.GetX(),
                y=atom.GetY(),
                z=atom.GetZ()
            ))

        return new_molecule
    except Exception as e:
        if node_runner:
            node_runner.error(f"Error reading CDX file {file_path}: {e}")
        raise ValueError(f"Invalid CDX file format or conversion error: {e}")

async def main():
    # Example usage (requires a real .cdx file to run)
    await context.initialize()
    cdx_path = Path().cwd().parent.parent / "tests" / "data" / "benzene.cdx"
    test_file = FileStack.from_local_file(cdx_path)
    mol = cdx_to_molecule(test_file)
    print(mol)


if __name__ == "__main__":
    asyncio.run(main())
