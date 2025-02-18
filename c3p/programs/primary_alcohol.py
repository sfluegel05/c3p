"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is defined as a hydroxyl group (-OH) that is attached to a saturated carbon atom, 
    which is bound to either three hydrogen atoms or one other carbon atom and two hydrogen atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise.
        str: Reason for classification.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define primary alcohol pattern (CH2OH or CH3OH)
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4H2][OX2H]")  # CH2OH where C is saturated
    simple_methanol_pattern = Chem.MolFromSmarts("[CH3][OX2H]")     # CH3OH
    
    # Check for the primary alcohol or methanol structure
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Contains a primary alcohol group"
    
    if mol.HasSubstructMatch(simple_methanol_pattern):
        return True, "Contains a methanol-like primary alcohol group"
    
    return False, "Does not contain a primary alcohol group"