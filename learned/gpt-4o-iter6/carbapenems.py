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
    
    # Define a refined SMARTS pattern for the carbapenem core:
    # This pattern captures a [3.2.0] bicyclic with a 4-membered beta-lactam ring
    # Allowing for substitutions at key positions
    core_pattern = Chem.MolFromSmarts("C1C2C=CC1C(=O)N2")  # Bicyclic structure for carbapenems

    # Additional patterns capturing common carbapenem features
    additional_patterns = [
        Chem.MolFromSmarts("C(=O)O"),  # Carboxyl group typically present in carbapenems
        Chem.MolFromSmarts("CS"),      # Thioether group, common in side chains
        Chem.MolFromSmarts("CC(=O)N"), # Amide linkage
    ]

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the core carbapenem pattern in the molecule
    if not mol.HasSubstructMatch(core_pattern):
        return False, "Does not match the carbapenem core structure"

    # Check for additional common substructures
    for pattern in additional_patterns:
        if not mol.HasSubstructMatch(pattern):
            return False, "Missing common functional groups for carbapenems"

    # If core and additional features are present, classify as carbapenem
    return True, "Contains characteristic carbapenem core and common functional groups"