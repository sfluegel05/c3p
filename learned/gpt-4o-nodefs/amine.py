"""
Classifies: CHEBI:32952 amine
"""
from rdkit import Chem

def is_amine(smiles: str):
    """
    Determines if a molecule is an amine based on its SMILES string.
    An amine contains a nitrogen atom with at least one alkyl or aryl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define amine SMARTS patterns
    primary_amine_pattern = Chem.MolFromSmarts("[NH2]")
    secondary_amine_pattern = Chem.MolFromSmarts("[NH]([C,R])")
    tertiary_amine_pattern = Chem.MolFromSmarts("[N]([C,R])([C,R])")
    
    # Check for amine groups
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains primary amine group (NH2)"
    if mol.HasSubstructMatch(secondary_amine_pattern):
        return True, "Contains secondary amine group (NH)"
    if mol.HasSubstructMatch(tertiary_amine_pattern):
        return True, "Contains tertiary amine group (N)"
    
    return False, "No amine functional group found"