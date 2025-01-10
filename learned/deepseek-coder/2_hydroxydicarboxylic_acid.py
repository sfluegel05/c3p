"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid
Definition: Any dicarboxylic acid carrying a hydroxy group on the carbon atom at position alpha to the carboxy group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for exactly two carboxylic acid groups (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) != 2:
        return False, f"Found {len(carboxylic_acid_matches)} carboxylic acid groups, need exactly 2"

    # Look for a hydroxyl group (-OH) on the alpha carbon (carbon adjacent to the carboxylic acid carbon)
    alpha_hydroxyl_pattern = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[OX2H1]")
    alpha_hydroxyl_matches = mol.GetSubstructMatches(alpha_hydroxyl_pattern)
    if len(alpha_hydroxyl_matches) == 0:
        return False, "No hydroxyl group found on the alpha carbon"

    # Ensure the hydroxyl group is on the alpha carbon of one of the carboxylic acid groups
    # This pattern ensures that the hydroxyl group is on the alpha carbon of one of the carboxylic acid groups
    specific_pattern = Chem.MolFromSmarts("[CX4][CX3](=[OX1])[OX2H1].[CX4][CX3](=[OX1])[OX2H1]")
    specific_matches = mol.GetSubstructMatches(specific_pattern)
    if len(specific_matches) == 0:
        return False, "Hydroxyl group not found on the alpha carbon of a carboxylic acid group"

    return True, "Contains two carboxylic acid groups and a hydroxyl group on the alpha carbon of one of them"