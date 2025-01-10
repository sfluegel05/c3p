"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid is a steroid with a ketone at position 3 and alpha configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid backbone (17 carbons in a fused ring system)
    steroid_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2[C@H]3[C@H]4CC[C@H]5CC(=O)CC[C@]5(C)[C@H]4CC[C@]3(C)[C@H]2CC1")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3-oxo group (ketone at position 3)
    oxo_pattern = Chem.MolFromSmarts("[CX3](=O)")
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    if not any(match[0] == 2 for match in oxo_matches):  # Assuming position 3 is index 2
        return False, "No 3-oxo group found"

    # Check for alpha configuration at position 5
    alpha_5_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2[C@H]3[C@H]4CC[C@H]5CC(=O)CC[C@]5(C)[C@H]4CC[C@]3(C)[C@H]2CC1")
    if not mol.HasSubstructMatch(alpha_5_pattern):
        return False, "No alpha configuration at position 5"

    return True, "Contains steroid backbone with 3-oxo group and alpha configuration at position 5"