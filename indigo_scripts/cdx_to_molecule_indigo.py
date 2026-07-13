from pathlib import Path
import logging

from molecular_qm_models import Atom, Molecule
from simstack.core.node import node
from simstack.models import Parameters, BooleanData
from simstack.models.files import FileStack

logger = logging.getLogger("CDX_to_Molecule")

@node(parameters=Parameters(resource="int-nano"))
def cdx_to_molecule_indigo(local_file_stack: FileStack, add_hydrogen: BooleanData, **kwargs) -> Molecule:
    """
    Convert a ChemDraw CDX file to a Molecule object using Indigo.
    
    :param local_file_stack: FileStack for the CDX file
    :param add_hydrogen: Whether to add hydrogen atoms to the molecule (default: True)
    :return: Molecule object
    """
    from indigo import Indigo
    task_id = kwargs.get("task_id", "NA")
    logger.info("Starting cdx-to-indigo V0.1")
    try:
        local_file = local_file_stack.get(Path(".."))
        indigo = Indigo()
        
        # Indigo can load CDX files directly
        mol = indigo.loadMoleculeFromFile(str(local_file))
        
        if add_hydrogen:
            mol.addImplicitHydrogens()
            
        # Convert Indigo molecule to SimStack Molecule
        simstack_mol = Molecule()
        for atom in mol.iterateAtoms():
            pos = atom.xyz()
            simstack_mol.add_atom(Atom(element=atom.symbol(), x=pos[0], y=pos[1], z=pos[2]))
            
        return simstack_mol

    except Exception as e:
        logger.error(f"Error converting CDX file to molecule with Indigo: {e}")
        raise ValueError(f"Failed to convert CDX file with Indigo: {e}")
