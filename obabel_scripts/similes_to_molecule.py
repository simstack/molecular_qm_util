import asyncio
import shutil
import subprocess
import tempfile
from pathlib import Path

from simstack.models import StringData, Parameters
from simstack.core.node import node
from molecular_qm_models import Molecule, Atom
from simstack.core.context import context

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


@node(parameters=Parameters(resource="int-nano"),force_rerun=True)
def smiles_to_molecule(smiles: StringData,**kwargs) -> Molecule:
    """
    Converts a SMILES string into a 3D Molecule object by invoking Open Babel's
    `obabel` as an external command (no Python OpenBabel bindings required).
    """
    node_runner = kwargs.get("node_runner", None)
    node_runner.info(f"Converting SMILES '{smiles}' to 3D coordinates")
    if isinstance(smiles, StringData):
        smiles = smiles.value
    if not smiles or not str(smiles).strip():
        raise ValueError("SMILES string cannot be empty")

    smiles = str(smiles).strip()

    #obabel = WindowsPath("C:/Program Files/OpenBabel-3.1.1") / "obabel.exe" # shutil.which("obabel")
    obabel = shutil.which("obabel")
    if not obabel:
        raise RuntimeError(
            "Could not find 'obabel' on PATH. Install Open Babel and ensure the 'obabel' "
            "executable is available on the command line."
        )

    # Use temp files to avoid shell escaping issues and to keep Windows happy.
    with tempfile.TemporaryDirectory(prefix="obabel_") as tmpdir:
        tmpdir_path = Path(tmpdir)
        smi_path = tmpdir_path / "input.smi"
        xyz_path = tmpdir_path / "output.xyz"

        smi_path.write_text(smiles + "\n", encoding="utf-8")

        # Notes:
        # - -ismi reads SMILES
        # - -oxyz writes XYZ
        # - --gen3d generates 3D coordinates
        # - -h adds hydrogens
        # - --minimize attempts a quick forcefield optimization
        # - --ff MMFF94 prefers MMFF94 (common); if not available, obabel may fall back or fail
        cmd = [
            obabel,
            "-ismi",
            str(smi_path),
            "-oxyz",
            "--gen3d",
            "-h",
            "--minimize",
            "--ff",
            "MMFF94",
            "--steps",
            "50",
            "-O",
            str(xyz_path),
        ]

        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            check=False,
        )

        if proc.returncode != 0:
            stderr = (proc.stderr or "").strip()
            stdout = (proc.stdout or "").strip()
            details = stderr if stderr else stdout
            raise RuntimeError(
                f"Open Babel (obabel) failed to convert SMILES '{smiles}'. "
                f"Exit code: {proc.returncode}. Details: {details}"
            )

        if not xyz_path.exists() or xyz_path.stat().st_size == 0:
            stderr = (proc.stderr or "").strip()
            stdout = (proc.stdout or "").strip()
            raise RuntimeError(
                f"Open Babel (obabel) produced no XYZ output for SMILES '{smiles}'. "
                f"stdout: {stdout} stderr: {stderr}"
            )

        xyz_text = xyz_path.read_text(encoding="utf-8", errors="replace")
        return _parse_xyz(xyz_text, smiles=smiles)

async def main():
    await context.initialize()
    smiles = StringData(value="C1=CC=CC=C1")
    mol = smiles_to_molecule(smiles)
    print(mol)

if __name__ == "__main__":
    asyncio.run(main())