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

    # Revised steroid pattern for four ring structure, allowing for varying unsaturations
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3=CC=CC4')  # Improved representation
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Revised 3beta-hydroxy group pattern with emphasis on specific stereochemistry
    # Ensure pattern covers the beta orientation on carbon 3 of the steroid ring
    hydroxy_beta_pattern = Chem.MolFromSmarts('[C@@H]1(C)O[*:2]')  # Beta hydroxyl at C3
    if not mol.HasSubstructMatch(hydroxy_beta_pattern):
        return False, "3beta-hydroxy group not properly oriented or absent"

    return True, "Molecule correctly classified as a 3beta-hydroxy steroid"