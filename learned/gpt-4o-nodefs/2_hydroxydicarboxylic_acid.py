"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid typically has a hydroxyl group on the second carbon
    and includes two carboxylic acid groups along the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create a SMARTS pattern for the 2-hydroxy group.
    hydroxy_pattern = Chem.MolFromSmarts("[#6][C@H](O)[#6]")
    
    # Create a SMARTS pattern for two separate carboxylic acid groups with some flexibility.
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX1H]")

    # Look for the 2-hydroxy group (OH on the second carbon or equivalent pattern).
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 2-hydroxy group found at an expected position"

    # Match two separate carboxylic acid groups
    carboxyl_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if len(carboxyl_matches) != 2:
        return False, f"Found {len(carboxyl_matches)} carboxylic acid groups, need exactly 2"

    # Apply further checks on the distance and separation logic if the logic is too prone to matches.
    # This can include bond count checks or heavier use of the compound connectivity map.
    # Currently, as per initial understanding focusing on the substructure match approach.
    
    return True, "Contains 2-hydroxy group on second position and two separate carboxylic acid groups"