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

    # Detect a spiro junction: two rings sharing one atom
    spiro_pattern = Chem.MolFromSmarts("[*]@[*]")
    if not mol.HasSubstructMatch(spiro_pattern):
        return False, "No spiro junction found"

    # Detect cyclic ketal: O-C-O in a ring context
    ketal_pattern = Chem.MolFromSmarts("[O][C]([O])([O])")
    if not mol.HasSubstructMatch(ketal_pattern):
        return False, "No cyclic ketal found"
    
    return True, "Contains a spiro junction with a cyclic ketal group"