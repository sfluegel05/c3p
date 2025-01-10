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

    # Look for steroid tetracyclic backbone: Allow for more flexibility with single/double bonds
    steroid_pattern = Chem.MolFromSmarts('C1C2C3C4CCCC(C4)C3C(C2)C1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Look for 3beta-hydroxy group - ensure beta configuration (stereochemistry specifics)
    # We assume correct sterochemical representations in the input SMILES for beta positioning
    hydroxy_beta_pattern = Chem.MolFromSmarts('C[C@H](O)[C@@H]')
    if not mol.HasSubstructMatch(hydroxy_beta_pattern):
        return False, "3beta-hydroxy group not found or not in beta orientation"

    return True, "Molecule classified as a 3beta-hydroxy steroid"

# This function now aims to more effectively identify core steroid structure
# while accounting for the orientation of the hydroxyl group using refined SMARTS patterns.