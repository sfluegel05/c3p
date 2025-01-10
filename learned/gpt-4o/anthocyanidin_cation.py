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
    
    # Check for the flavylium cation core structure which includes the [o+] and pyran ring
    # This is a more general pattern for aromatic chromenylium structure with [o+] core
    flavylium_pattern = Chem.MolFromSmarts('c1c[o+]c2cc(O)c(O)cc2c(O)c1') # Tentative SMARTS for flavylium
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium cation core structure found"
    
    # Identify the presence of multiple hydroxyl (OH) or methoxy (OMe) groups on the flavonoid rings
    hydroxyl_methoxy_pattern = Chem.MolFromSmarts('[cH]-[OH,OMe]')
    hydroxyl_methoxy_count = len(mol.GetSubstructMatches(hydroxyl_methoxy_pattern))
    if hydroxyl_methoxy_count < 3:  # Assume 3 such groups for more specificity in anthocyanidins
        return False, f"Found {hydroxyl_methoxy_count} hydroxyl/methoxy groups, require at least 3"

    return True, "Contains characteristic flavylium cation core structure with sufficient hydroxyl/methoxy groups"