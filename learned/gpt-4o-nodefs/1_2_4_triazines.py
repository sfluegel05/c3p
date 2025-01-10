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

    # General 1,2,4-triazine pattern allowing substituents
    triazine_patterns = [
        Chem.MolFromSmarts("n1cncnc1"),  # Basic 1,2,4-triazine ring with variable substituents
        Chem.MolFromSmarts("n1c[nH]c[nH]c1")  # another potential protonated variant, etc.
        # Add more patterns if needed for other structural variants
    ]

    for pattern in triazine_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a 1,2,4-triazine ring structure"

    return False, "Does not contain a 1,2,4-triazine ring structure"