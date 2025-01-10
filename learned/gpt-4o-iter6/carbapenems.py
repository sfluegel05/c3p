"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determine if a given SMILES string is a carbapenem.

    Args:
        smiles (str): The SMILES string of the chemical entity.
    
    Returns:
        bool: True if SMILES represents a carbapenem, False otherwise.
        str: Explanation of the reasoning for the classification.
    """
    
    # Define the core carbapenem pattern:
    # Bicyclic system with the beta-lactam ring
    core_pattern = Chem.MolFromSmarts("C1C(C)C2C1N(C2=O)C(=O)C")  # Captures [3.2.0] with beta-lactam

    # Additional patterns capturing variations in substitution that retain carbapenem identity
    additional_patterns = [
        Chem.MolFromSmarts("C(=O)O"),  # Carboxylic acid group
        Chem.MolFromSmarts("[S]"),  # Presence of a sulfur atom, often in thioether form
    ]

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Check for the core carbapenem structure
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match the carbapenem core structure"

    # Check for characteristic functional groups or atoms
    for pattern in additional_patterns:
        if not mol.HasSubstructMatch(pattern):
            return False, "Missing common functional group"

    # If the characteristic features are present, classify as carbapenem
    return True, "Contains characteristic carbapenem core and functional groups"