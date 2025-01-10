"""
Classifies: CHEBI:27300 vitamin D
"""
from rdkit import Chem

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for secosteroid structure: look for a broken B-ring pattern
    # General simplified pattern for secosteroids could be further improved
    seco_steroid_pattern = Chem.MolFromSmarts("C1CCC2=C(CC1)CCC3C2CCC4[C@@H]3[C@@H](CCCC4)C")
    if not mol.HasSubstructMatch(seco_steroid_pattern):
        return False, "No secosteroid core structure found"
        
    # Check for hydroxyl groups at specific vitamin D positions (3 position is critical)
    # SMARTS for hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[CX4][OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if len(hydroxyl_matches) < 2:  # Vitamin D3 generally has two primary hydroxyl groups
        return False, "Not enough hydroxyl groups, found less than 2"
    
    # Check for the presence of a conjugated triene system, specific to the vitamin D3 structure
    triene_pattern = Chem.MolFromSmarts("C=C/C=C/C=C")
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene system detected for Vitamin D"

    return True, "Matches vitamin D core structure and functional groups"