"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is defined as a hydroxyl group (-OH) attached to a saturated carbon,
    which is bonded to either two hydrogen atoms and one carbon atom or three hydrogen atoms.
    
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
    
    # Define a primary alcohol pattern considering different cases
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4;H2][CH2][OX2H]")  # Primary alcohol where carbon has two H
    methanol_pattern = Chem.MolFromSmarts("[CH3][OX2H]")                # Methanol (special case)

    # Check for primary alcohol substructure
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Contains a primary alcohol group"
    
    # Check for methanol-like structure
    if mol.HasSubstructMatch(methanol_pattern):
        return True, "Contains a methanol-like primary alcohol group"
    
    return False, "Does not contain a primary alcohol group"