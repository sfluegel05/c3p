"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is a six-membered heterocyclic ring containing three nitrogens at
    positions 1, 2, and 4. This pattern may exist with various substituents attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1,2,4-triazine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS patterns for identifying 1,2,4-triazine ring
    triazine_patterns = [
        Chem.MolFromSmarts("n1ncncc1"),  # 1,2,4-triazine minimal pattern
        Chem.MolFromSmarts("c1ncncn1"),  # aromatic representation
        Chem.MolFromSmarts("n1c[nH]cnc1"),  # potential tautomer/protonated form
        Chem.MolFromSmarts("n1[nH]cc[nH]c1")  # another possible variant
    ]

    # Check for a match with any of the patterns
    for pattern in triazine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a 1,2,4-triazine ring structure"

    return False, "Does not contain a 1,2,4-triazine ring structure"