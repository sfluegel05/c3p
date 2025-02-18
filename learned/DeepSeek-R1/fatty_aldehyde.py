"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: CHEBI:174450 fatty aldehyde
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde derived from a fatty acid, with a carbonyl group at one end.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for exactly one aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if len(aldehyde_matches) != 1:
        return False, f"Found {len(aldehyde_matches)} aldehyde groups (needs exactly 1)"

    # Check that the aldehyde carbon is terminal (only one adjacent carbon)
    aldehyde_carbon = aldehyde_matches[0][0]
    neighbors = [atom for atom in mol.GetAtomWithIdx(aldehyde_carbon).GetNeighbors() if atom.GetAtomicNum() == 6]
    if len(neighbors) != 1:
        return False, "Aldehyde group is not at the end of a carbon chain"

    # Check for absence of carboxylic acid groups
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Contains carboxylic acid group"

    return True, "Terminal aldehyde group with a carbon chain, derived from fatty acid"