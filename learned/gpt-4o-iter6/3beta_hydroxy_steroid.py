"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.
    A 3beta-hydroxy steroid is a steroid with a hydroxyl group in the beta position at carbon 3.

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
    
    # Define a more flexible steroid skeletal pattern
    steroid_pattern = Chem.MolFromSmarts('[#6]1[#6][#6]2[#6]([#6]1)[#6][#6]3[#6]([#6]2)[#6][#6]4[#6]=[#6][#6][#6]4[#6]3')  # Capturing more variability in structure
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3beta-hydroxy group
    hydroxy_beta_pattern = Chem.MolFromSmarts('[C@H](O)[C@@H]1CC[C@H]2[C@@H]1CCC3[C@@H]2[C@H](C1=CC=[C@H]4CC[C@H](O)C=C41)[C@@H]3')  # Ensure correct stereochemistry for 3beta-hydroxy
    if not mol.HasSubstructMatch(hydroxy_beta_pattern):
        return False, "3beta-hydroxy group not properly oriented or absent"

    return True, "Molecule correctly classified as a 3beta-hydroxy steroid"