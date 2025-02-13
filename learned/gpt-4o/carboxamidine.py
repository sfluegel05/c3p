"""
Classifies: CHEBI:35359 carboxamidine
"""
from rdkit import Chem

def is_carboxamidine(smiles: str):
    """
    Determines if a molecule is a carboxamidine based on its SMILES string.
    A carboxamidine is characterized by the presence of a -C(=N)-NR2 group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carboxamidine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxamidine SMARTS pattern
    carboxamidine_pattern = Chem.MolFromSmarts("C(=N)[N&R]")  # Simplified pattern for carboxamidine
    if mol.HasSubstructMatch(carboxamidine_pattern):
        return True, "Contains carboxamidine group (-C(=N)-NR2)"

    return False, "Does not contain a carboxamidine group"

# Example usage:
# test_smiles = "NC(=N)c1ccccc1"  # benzamidine
# is_carboxamidine(test_smiles)