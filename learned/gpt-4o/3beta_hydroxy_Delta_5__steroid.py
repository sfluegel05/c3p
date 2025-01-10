"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid contains a 3beta-hydroxy group and a double bond between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a 3beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C@@H]1CC[C@@H](C2=CC=C3[C@H](C=C2[C@H]1)C(C3)=O)")
    
    # Define SMARTS pattern for a Delta(5) double bond (C=C) within a steroid structure
    delta5_pattern = Chem.MolFromSmarts("[C:5]=[C:6]")
    
    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group found"
        
    # Check for Delta(5) double bond
    if not mol.HasSubstructMatch(delta5_pattern):
        return False, "No Delta(5) double bond found (C5-C6 double bond)"

    return True, "Contains 3beta-hydroxy group with Delta(5) double bond (C5-C6 double bond)"