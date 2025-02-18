"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is defined as a hydroxyl group (-OH) that is attached to a carbon 
    which is only bound to either two hydrogen atoms and one carbon atom or three hydrogen atoms,
    or considered the methanol form of CH3OH.
    
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
    
    # Define comprehensive primary alcohol pattern
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4;!$(C([#6])[#6])][CH2][OX2H]")  # Keep the CH2 not bonded to 2 other carbons 
    methanol_pattern = Chem.MolFromSmarts("[CH3][OX2H]")  # Explicitly check methanol

    # Check if the molecule matches either pattern
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Contains a primary alcohol group"
    
    if mol.HasSubstructMatch(methanol_pattern):
        return True, "Contains a methanol-like primary alcohol group"
    
    return False, "Does not contain a primary alcohol group"