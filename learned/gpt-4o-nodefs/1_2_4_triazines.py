"""
Classifies: CHEBI:39410 1,2,4-triazines
"""
from rdkit import Chem

def is_1_2_4_triazines(smiles: str):
    """
    Determines if a molecule is a 1,2,4-triazine based on its SMILES string.
    A 1,2,4-triazine is a six-membered heterocyclic ring with three nitrogens at 
    positions 1, 2, and 4.

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

    # SMARTS patterns for identifying 1,2,4-triazine ring structures
    triazine_patterns = [
        Chem.MolFromSmarts("[nR]1cncnc1"),  # Typical 1,2,4-triazine pattern
        Chem.MolFromSmarts("c1ncncn1"),    # Aromatic 1,2,4-triazine
        Chem.MolFromSmarts("n1cncnc1"),    # Non-aromatic variant
        Chem.MolFromSmarts("c1ncn[nH]c1"), # Potential protonated/triazinium form
        Chem.MolFromSmarts("n1c[nH]cnc1")  # Potential tautomer
    ]

    # Check for a match with any of the patterns
    for pattern in triazine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a 1,2,4-triazine ring structure"

    return False, "Does not contain a 1,2,4-triazine ring structure"