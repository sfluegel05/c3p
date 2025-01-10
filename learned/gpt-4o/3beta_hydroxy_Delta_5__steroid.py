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
    
    # SMARTS pattern for a 3beta-hydroxy group on a steroid backbone
    hydroxy_pattern = Chem.MolFromSmarts("[C@H](O)[C@@H]1CC[C@H](C)C2=CC=C3[C@@H](CCCC3)C[C@H]12") # Example steroid backbone with 3β-hydroxy
    
    # SMARTS pattern for a Delta(5) double bond within a steroid structure
    delta5_pattern = Chem.MolFromSmarts("C1=CCC2C1CC[C@@H]3CCC=C4C=C[C@H]4[C@H]23") # Correct position for Delta(5) within a typical steroid
    
    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group found"

    # Check for Delta(5) double bond
    if not mol.HasSubstructMatch(delta5_pattern):
        return False, "No Delta(5) double bond found (C5-C6 double bond)"

    return True, "Contains 3beta-hydroxy group with Delta(5) double bond (C5-C6 double bond)"