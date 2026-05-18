import subprocess
from pathlib import Path
from molecular_qm_models import Molecule

def smiles_to_molecule(smiles_string: str) -> Molecule:
    # uv run will create a cached, isolated env just for this call
    script_path = Path(__file__).parent /  "rd_kit_scripts" / "generate_3d_structure_from_smiles.py"
    result = subprocess.run(
        ["uv", "run", script_path, smiles_string],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise RuntimeError(f"Error running script: {result.stderr}")

    new_molecule = Molecule.from_sdf(result.stdout)
    return new_molecule