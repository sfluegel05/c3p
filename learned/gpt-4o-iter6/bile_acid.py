"""
Classifies: CHEBI:3098 bile acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    Bile acids are hydroxy-5beta-cholanic acids that have a steroid nucleus and specific functional groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General pattern for 5beta-cholanic acid framework (allowing flexibility in peripherals)
    sterane_pattern = Chem.MolFromSmarts('[C@]12[C@@H]3CC[C@]4([C@H]3C=C[C@]4([C@@H]2CC[C@H]1[CH3])C)CC(=O)')
    if not mol.HasSubstructMatch(sterane_pattern):
        return False, "5-beta-steroid nucleus not identified"

    # Check for hydroxyl groups
    hydroxy_pattern = Chem.MolFromSmarts('C(O)')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxyl_matches) < 1:
        return False, "Insufficient number of hydroxy groups (less than one)"

    # Check for the carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)[O-]')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Carboxylic acid group not identified"

    return True, "Valid bile acid: contains flexible steroid framework, hydroxyl groups, and a terminal carboxylic group"

# Example usage:
# result, reason = is_bile_acid("O[C@@H]1[C@]2([C@]...")