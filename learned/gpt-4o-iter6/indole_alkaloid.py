"""
Classifies: CHEBI:38958 indole alkaloid
"""
from rdkit import Chem

def is_indole_alkaloid(smiles: str):
    """
    Determines if a molecule is an indole alkaloid based on its SMILES string.
    An indole alkaloid must have an indole skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an indole alkaloid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for indole core
    indole_pattern = Chem.MolFromSmarts('C1=CC2=C(C=C1)NC=C2')
    if not mol.HasSubstructMatch(indole_pattern):
        return False, "No indole skeleton found"
    
    # Additional checks for alkaloid properties
    # Often containing extra nitrogen or specific functional groups
    
    # Define a basic nitrogen containing group pattern
    nitrogen_pattern = Chem.MolFromSmarts('[#7]')
    if not mol.HasSubstructMatch(nitrogen_pattern):
        return False, "Lacks additional nitrogen typically found in alkaloids"
    
    # If more specific functional groups need to be checked, define here

    return True, "Contains indole skeleton and additional features typical of alkaloids"

# Usage example for ergoptine
result, reason = is_indole_alkaloid("O=C1N2[C@H]([C@@]3(O[C@](C(N3[C@H]1CC(C)C)=O)(NC(=O)[C@@H]4C=C5C6=C7C(NC=C7C[C@H]5N(C4)C)=CC=C6)CC)O)CCC2")
print(result, reason)