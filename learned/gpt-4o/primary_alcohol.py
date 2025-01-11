"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is characterized by a hydroxyl (-OH) group attached
    to a carbon which is attached to at least two hydrogen atoms, or a
    carbon with one other carbon and two hydrogen atoms.

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
    
    # Define an improved SMARTS pattern that captures a broader definition of primary alcohol
    # The main features of a primary alcohol: OH group not linked to carbon with more than one non-hydrogen substituent
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4;H2,H1][OH]")

    # Check the pattern match
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "SMILES indicates a primary alcohol structure"
    else:
        return False, "No primary alcohol structure found in SMILES"

# Example of usage:
# smiles = '...'  # Replace with actual SMILES string
# is_primary_alcohol(smiles)