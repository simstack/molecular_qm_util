import logging
import pubchempy as pcp
from collections import Counter
from molecular_qm_models import Molecule

logger = logging.getLogger(__name__)

def compute_iupac_name(molecule: "Molecule") -> str:
    """
    Compute the IUPAC name for a molecule using PubChemPy.
    Falls back to molecular formula if naming fails.

    :param molecule: A Molecule object.
    :return: IUPAC name string or molecular formula.
    """
    if not molecule.atoms:
        raise ValueError("Cannot compute IUPAC name for empty molecule")

    try:
        # Get SMILES representation to query PubChem
        smiles = molecule.smiles or molecule.make_smiles()
        
        # Query PubChem
        compounds = pcp.get_compounds(smiles, namespace='smiles')
        
        if compounds and compounds[0].iupac_name:
            return compounds[0].iupac_name
        else:
            logger.warning(f"No IUPAC name found on PubChem for SMILES: {smiles}")

    except Exception as e:
        logger.error(f"Error retrieving IUPAC name from PubChem: {e}")

    # Fallback to alphabetically ordered formula
    element_counts = Counter(atom.element for atom in molecule.atoms)
    sorted_elements = sorted(element_counts.keys())
    return ''.join(f"{elem}{element_counts[elem]}" for elem in sorted_elements)
