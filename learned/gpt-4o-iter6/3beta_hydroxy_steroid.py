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

    # Improved steroid tetracyclic backbone pattern
    # Accounting for variability in unsaturation, stereochemistry
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4(C3=CC=CC4)')  # Flexible circular and linear cases
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Look for 3beta-hydroxy group - ensure beta configuration (using chiral specific pattern)
    # Utilizing more flexible chiral specifications for hydroxyl at C3
    hydroxy_beta_pattern = Chem.MolFromSmarts('C[C@@H](O)[C]')  # Beta hydroxyl group needed
    if not mol.HasSubstructMatch(hydroxy_beta_pattern):
        return False, "3beta-hydroxy group not properly oriented or absent"
    
    return True, "Molecule correctly classified as a 3beta-hydroxy steroid"

# This code aims to more accurately match inherent structural diversity present within these steroids while 
# focusing on precise detection of functional groups and configurational isomerism.