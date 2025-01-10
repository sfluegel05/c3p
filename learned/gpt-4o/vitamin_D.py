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

    # Broaden the search for secosteroid structure, looking for typical seco-B-ring pattern
    # This pattern is indicative but flexible for ring openings
    open_b_ring_pattern = Chem.MolFromSmarts("C1CCC2C(C1)CC=C3C2CCCC3")
    if not mol.HasSubstructMatch(open_b_ring_pattern):
        return False, "No secosteroid-like open B-ring structure found indicative of vitamin D"

    # Check for required triene system, considering extended conjugation possibilities
    triene_conjugated_pattern = Chem.MolFromSmarts("C=C/C=C/C=C")
    if not mol.HasSubstructMatch(triene_conjugated_pattern):
        return False, "No conjugated triene system detected, typical of vitamin D"

    # Check for required hydroxyl groups presence
    # Vitamin D compounds should have variations of these groups, allow any position with leniency here
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    if len(hydroxyl_matches) < 1:  # Allow at least one hydroxyl group (original or derivatized)
        return False, "No hydroxyl group(s) detected, typical of vitamin D"

    return True, "Matches typical vitamin D structural features"