"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are characterized by a positively-charged flavylium core with hydroxyl or methoxy groups.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is an anthocyanidin cation, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the flavylium cation core structure which includes the [o+] charge and chromenylium structure
    # Using a broader SMARTS pattern to capture various flavylium configurations with appropriate oxygen charge
    flavylium_pattern = Chem.MolFromSmarts('c1cc2[o+]ccc2cc1') # A more generic representation
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium cation core structure found"
    
    # Identify and count the presence of hydroxyl (OH) or methoxy (OMe) groups
    hydroxyl_methoxy_pattern = Chem.MolFromSmarts('[OH,OMe]')
    hydroxyl_methoxy_count = len(mol.GetSubstructMatches(hydroxyl_methoxy_pattern))
    if hydroxyl_methoxy_count < 3:  # Expecting at least 3 such groups for anthocyanidin characterization
        return False, f"Found {hydroxyl_methoxy_count} hydroxyl/methoxy groups, require at least 3"

    return True, "Contains characteristic flavylium cation core structure with sufficient hydroxyl/methoxy groups"