"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is defined by a six-membered aromatic ring with nitrogen
    atoms at positions 1, 2, and 4 of the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more flexible SMARTS patterns for 1,2,4-triazine
    triazine_patterns = [
        Chem.MolFromSmarts("n1cnnc[nH]1"),  # without aromaticity
        Chem.MolFromSmarts("n1[nH]cnn1"),  # oxidized form
        Chem.MolFromSmarts("c1nncnc1"),    # full aromatic form
        Chem.MolFromSmarts("n1cnncn1")     # original pattern
    ]

    # Check for each 1,2,4-triazine pattern in molecule
    for pattern in triazine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains 1,2,4-triazine ring structure"

    # If none of the patterns match
    return False, "Does not contain 1,2,4-triazine ring structure"