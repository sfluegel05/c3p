"""
Classifies: CHEBI:15734 primary alcohol
"""
"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol has a hydroxy group (-OH) attached to a carbon atom that is either:
    - bonded to three hydrogen atoms (CH3-OH), or
    - bonded to one other carbon atom and two hydrogen atoms (R-CH2-OH).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure pattern for a primary alcohol
    # Pattern 1: -CH2-OH (R-CH2-OH)
    pattern1 = Chem.MolFromSmarts("[CH2][OH]")
    # Pattern 2: -CH3-OH (CH3-OH)
    pattern2 = Chem.MolFromSmarts("[CH3][OH]")

    # Check if the molecule matches either pattern
    if mol.HasSubstructMatch(pattern1):
        return True, "Contains a -CH2-OH group (primary alcohol)"
    elif mol.HasSubstructMatch(pattern2):
        return True, "Contains a -CH3-OH group (primary alcohol)"
    else:
        return False, "No primary alcohol group (-CH2-OH or -CH3-OH) found"