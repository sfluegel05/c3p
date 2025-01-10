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
    
    # Define a more precise SMARTS pattern for the carbapenem core:
    # This pattern involves a [3.2.0] bicyclic system with a 4-membered beta-lactam ring
    core_pattern = Chem.MolFromSmarts("[CX4]1=[NX3][C@@H]2C[C@H]1C(=O)N2")  # Adjusted pattern

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the core carbapenem pattern in the molecule
    if mol.HasSubstructMatch(core_pattern):
        # Additional checks could be included to verify specific substitutions or functionalities if required
        return True, "Matches the carbapenem core structure"
    else:
        return False, "Does not match the carbapenem core structure"