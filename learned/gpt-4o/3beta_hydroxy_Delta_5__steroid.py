"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid is a 3beta-hydroxy-steroid that contains a double bond
    between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 3beta-hydroxy group pattern (OH attached to C-3)
    hydroxy_3beta_pattern = Chem.MolFromSmarts("[C@@H](O)[C@]")

    # Steroid nucleus pattern for fused ABCD rings
    steroid_pattern = Chem.MolFromSmarts("C1CC2CCC3C4CCC(C)(C=CC4=CC3=CC2=C1)")

    # Delta(5) double bond pattern at specific positions
    delta5_pattern = Chem.MolFromSmarts("C1=CC=CC2CCC3C4CCC(C=CC4=CC3=CC2=CC1)")

    # Check for steroid backbone
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected."

    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_3beta_pattern):
        return False, "No 3beta-hydroxy group match found."

    # Check for specific Delta(5) double bond
    if not mol.HasSubstructMatch(delta5_pattern):
        return False, "No specific Delta(5) double bond match found in ring."

    return True, "Molecule is classified as 3beta-hydroxy-Delta(5)-steroid."