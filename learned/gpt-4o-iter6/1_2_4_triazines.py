"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is characterized by nitrogen atoms at positions 1, 2, and 4
    of a six-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern allowing for aromatic and non-aromatic forms
    # Example pattern (broad, to capture nitrogens at positions): 
    triazine_patterns = [
        Chem.MolFromSmarts("n1cnnc1"),  # Non-aromatic core
        Chem.MolFromSmarts("c1ncnnc1"), # Aromatic core
        Chem.MolFromSmarts("n1c[nH]nc1"), # Mixed aromaticity core
        Chem.MolFromSmarts("c1nnc[nH]c1") # Alternative tautomeric/core form
    ]
    
    # Check if the molecule contains any of the 1,2,4-triazine substructures
    for pattern in triazine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains the 1,2,4-triazine core structure"
    
    return False, "Does not contain the 1,2,4-triazine core structure"