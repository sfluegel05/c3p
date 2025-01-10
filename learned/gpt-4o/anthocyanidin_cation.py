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
    
    # Check for flavylium core structure with positive charge ([o+] aromatic system)
    flavylium_cation_pattern = Chem.MolFromSmarts('c1cc2cc[c+]cc2[o]1')
    if not mol.HasSubstructMatch(flavylium_cation_pattern):
        return False, "No flavylium cation core structure found"
    
    # Identify presence of multiple hydroxyl (OH) groups attached to rings
    hydroxyl_pattern = Chem.MolFromSmarts('c[OH]')
    hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if hydroxyl_count < 2:
        return False, f"Found {hydroxyl_count} hydroxyl groups, require at least 2"

    return True, "Contains flavylium cation core structure with sufficient hydroxyl groups"