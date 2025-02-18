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

    # Improve primary alcohol pattern to support broader definitions
    primary_pattern = Chem.MolFromSmarts("[CH2][OX2H]")  # Primary alcohol with two hydrogens on C
    methanol_pattern = Chem.MolFromSmarts("[CH3][OX2H]")  # Methanol

    # Check for match using newly defined patterns
    if mol.HasSubstructMatch(primary_pattern):
        return True, "Contains a primary alcohol (CH2-OH) group"

    if mol.HasSubstructMatch(methanol_pattern):
        return True, "Contains a methanol (CH3-OH) group"

    return False, "Does not contain a primary alcohol group"