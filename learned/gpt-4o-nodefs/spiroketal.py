"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal should have a spiro junction and a cyclic ketal structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved detection of a spiro junction (an atom that is part of two rings)
    spiro_pattern = Chem.MolFromSmarts("[R]@[R]")  # Atoms that are part of two ring closures
    if not mol.HasSubstructMatch(spiro_pattern):
        return False, "No spiro junction found"

    # Improved detection of cyclic ketal (RCOCR')
    ketal_pattern = Chem.MolFromSmarts("O[C;R]O")  # Looking for C-O-C linkages within ring context
    if not mol.HasSubstructMatch(ketal_pattern):
        return False, "No cyclic ketal found"
    
    return True, "Contains a spiro junction with a cyclic ketal group"