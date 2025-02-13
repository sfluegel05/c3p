"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies oxo fatty acids based on their SMILES string.
"""
from rdkit import Chem

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    An oxo fatty acid is defined as any fatty acid containing at least one aldehydic or ketonic group
    in addition to the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group pattern (C(=O)O)
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for aldehydic or ketonic groups (C=O, distinct from the carboxyl group)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    
    # Check that there's more than one carbonyl group, confirming at least one as oxo group distinct from carboxyl
    if len(carbonyl_matches) <= 1:
        return False, f"Insufficient carbonyl groups, need at least one additional to carboxylic group"

    return True, "Contains carboxylic acid group and additional oxo group(s) (aldehydic/ketonic)"