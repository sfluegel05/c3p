"""
Classifies: CHEBI:16366 anthocyanidin cation
"""
from rdkit import Chem

def is_anthocyanidin_cation(smiles: str):
    """
    Determines if a molecule is an anthocyanidin cation based on its SMILES string.
    Anthocyanidin cations are aglycons of anthocyanins with a characteristic flavylium cation structure.
    
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
    
    # Check for the flavylium cation structure with aromatic pyran
    # General flavylium pattern can be complex due to possible substitutions
    flavylium_pattern = Chem.MolFromSmarts('c1cc(OC2=C(C=CC=C2)[O+]C=C1)c3ccccc3') # Possible flavylium cation core
    if not mol.HasSubstructMatch(flavylium_pattern):
        return False, "No flavylium cation core structure found"
    
    # Identify presence of hydroxyl (OH) or methoxy (OCH3) groups substituted onto the aromatic rings
    # Check for at least 2 such groups could help specify anthocyanidin
    hydroxyl_methoxy_pattern = Chem.MolFromSmarts('[cH]-[OH,OMe]')
    hydroxyl_methoxy_count = len(mol.GetSubstructMatches(hydroxyl_methoxy_pattern))
    if hydroxyl_methoxy_count < 2:
        return False, f"Found {hydroxyl_methoxy_count} hydroxyl/methoxy groups, require at least 2"

    return True, "Contains flavylium cation core structure with sufficient functional groups"