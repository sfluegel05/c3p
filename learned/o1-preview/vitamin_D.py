"""
Classifies: CHEBI:27300 vitamin D
"""
"""
Classifies: vitamin D
"""

from rdkit import Chem

def is_vitamin_D(smiles: str):
    """
    Determines if a molecule is a vitamin D compound based on its SMILES string.
    Vitamin D compounds are secosteroids with a broken B-ring (9,10-seco-steroids),
    characterized by a triene system in the opened ring and specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a vitamin D compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the vitamin D secosteroid skeleton SMARTS pattern
    # The pattern represents the opened B-ring with a conjugated triene system
    vitamin_d_skeleton_smarts = """
    C1=C[C@H]2CC[C@@]3(C1)C=C\C=C\C4=C3CCCC4
    """
    vitamin_d_pattern = Chem.MolFromSmarts(vitamin_d_skeleton_smarts.strip())
    if vitamin_d_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Check for the secosteroid skeleton match
    if not mol.HasSubstructMatch(vitamin_d_pattern):
        return False, "No vitamin D secosteroid skeleton found"
    
    # Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts('[OX2H]')
    hydroxyls = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyls) == 0:
        return False, "No hydroxyl groups found"
    
    # Additional check for secosteroid core with triene system
    triene_pattern = Chem.MolFromSmarts('C=C-C=C-C=C')
    if not mol.HasSubstructMatch(triene_pattern):
        return False, "No conjugated triene system found"

    return True, "Contains vitamin D secosteroid skeleton with triene system"