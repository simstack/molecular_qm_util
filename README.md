# molecular_qm_util
Utilities for molecular electronic structure calculations with Simstack

## Functions and their dependencies

The package is organized into sub-modules named after the third-party library
they rely on. The table below documents which package powers which function so
you can tell, at a glance, what needs to be installed for a given feature.

| Function | Import path | Backend package(s) | What it does |
| --- | --- | --- | --- |
| `compute_iupac_name` | `molecular_qm_util.pubchempy_scripts` | `pubchempy` | Looks up the IUPAC name of a `Molecule` on PubChem (falls back to a molecular formula when no name is found). |
| `smiles_to_molecule_obabel` (Open Babel) | `molecular_qm_util.obabel_scripts` | `openbabel` (`openbabel-wheel`), `simstack` | Converts a SMILES string to a 3D `Molecule` using Open Babel (adds hydrogens, builds 3D coordinates, MMFF94/UFF minimization). Runs as a SimStack node. |
| `cdx_to_molecule` | `molecular_qm_util.obabel_scripts` | `openbabel` (`openbabel-wheel`), `simstack` | Reads a ChemDraw `.cdx` file into a `Molecule` using Open Babel. |
| `molecule_to_smiles` / `compute_smiles` / `molecule_to_canonical_smiles` | `molecular_qm_util.rdkit_scripts` | `rdkit` | Converts a `Molecule` (3D coordinates) to a canonical SMILES string using RDKit bond perception. The three names are aliases of the same function. |
| `smiles_to_molecule_obabel` (RDKit) | `molecular_qm_util.rdkit_scripts` | `rdkit` | Converts a SMILES string to a 3D `Molecule` using RDKit ETKDG embedding. |
| `get_3d` | `molecular_qm_util.rdkit_scripts` | `rdkit` | Generates a 3D MolBlock from a SMILES string via RDKit ETKDG embedding. |
| `get_rotatable_bonds` | `molecular_qm_util.rdkit_scripts` | `rdkit` | Finds rotatable bonds (as dihedral atom quadruples) in an RDKit molecule. |
| `rotate_bond_rigidly` / `get_fragment_atoms` | `molecular_qm_util.rdkit_scripts` | `rdkit`, `numpy` | Rigidly rotates a bond to enumerate conformers around a dihedral. |

### Conversion helpers (`molecular_qm_util.util`)

| Function | Backend package(s) | What it does |
| --- | --- | --- |
| `rdkit_mol_to_molecule` | `rdkit`, `molecular_qm_models` | Converts an RDKit `Mol` (with a conformer) to a SimStack `Molecule`. |
| `pybel_mol_to_molecule` | `openbabel` (pybel), `molecular_qm_models` | Converts an Open Babel `pybel` molecule to a SimStack `Molecule`. |
| `simstack_molecule_to_rdkit` | `openbabel` (pybel), `rdkit`, `molecular_qm_models` | Converts a SimStack `Molecule` to an RDKit `Mol` via Open Babel. |

> Note: every function also uses `molecular_qm_models` for the `Molecule`/`Atom`
> data model. The `openbabel`- and `rdkit`-based helpers are imported lazily in
> `molecular_qm_util/__init__.py`, so a missing optional backend does not break
> importing `compute_iupac_name`.

