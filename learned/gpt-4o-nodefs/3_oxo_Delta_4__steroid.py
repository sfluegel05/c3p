"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule matches the class, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the steroid backbone (ABCD rings)
    steroid_pattern = Chem.MolFromSmarts("C1CCC2CCCC3C4CCC(=O)CC4CCC23")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Look for the 3-oxo group (C=O directly attached to a ring carbon)
    oxo_pattern = Chem.MolFromSmarts("C(=O)[C&R]")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "Missing 3-oxo group"
    
    # Look for the Delta(4) double bond (C=C between C4 and C5)
    delta4_pattern = Chem.MolFromSmarts("C=CC=O")
    if not mol.HasSubstructMatch(delta4_pattern):
        return False, "Missing Delta(4) double bond"

    return True, "Contains the structural features of a 3-oxo-Delta(4) steroid"

# Example of usage:
# result, reason = is_3_oxo_Delta_4__steroid("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1C[C@H](O)[C@@H]2O")
# print(result, reason)