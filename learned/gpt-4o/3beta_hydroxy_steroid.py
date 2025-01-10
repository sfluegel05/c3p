"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is characterized by a steroid structure with a hydroxyl
    group at the 3rd carbon in the beta orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a pattern for the steroid ABCD rings, allowing for variations in saturation
    steroid_pattern = Chem.MolFromSmarts("[#6]1-[#6]-[#6]-2-[#6]-[#6]3-[#6](-[#6]-[#6]-4-[#6R1]-[#6R1]([#6R1][#6]4)[#6R1]([#6R1][#6]3)[#6R1][#6]2)-[#6R1][#6]1")
    
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
    
    # Define a pattern for a 3beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@H](C)CC[C@]2(C)[C@@H]3CC[C@@H]4C([C@H](O)C)CCC3=[C@]2[C@H]4C1")
    
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group found"

    return True, "Contains steroid backbone with 3beta-hydroxy group"