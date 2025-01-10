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
    
    # Refined pattern for 5beta-cholanic acid framework
    cholanic_pattern = Chem.MolFromSmarts('[C@@H]1[C@H]([C@H]2[C@@H](CCC3[C@@H](C)C[C@H](O)[C@@H]4[C@@H]3CC[C@]4(C)[C@H]2CC1)O)C(C)C(=O)O')
    if not mol.HasSubstructMatch(cholanic_pattern):
        return False, "5beta-cholanic acid framework not found"
    
    # Check for at least 2 or 3 hydroxyl groups - this might depend on specific bile acid examples
    hydroxy_pattern = Chem.MolFromSmarts('[C@H](O)')
    hydroxyl_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Insufficient number of hydroxy groups (less than two)"

    # Check for the carboxylic acid group
    carboxylic_pattern = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "Carboxylic acid group not found"

    return True, "Valid bile acid: contains 5beta-cholanic framework, hydroxyl groups, and a terminal carboxylic group"

# Example usage:
# result, reason = is_bile_acid("O[C@@H]1[C@]2([C@]...")