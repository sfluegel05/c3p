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
    
    # Define a flexible steroid skeletal pattern with typical steroid rings (cyclic structure C17H28)
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3(O)CCC4') 
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Check for 3beta-hydroxy group, allowing for alternative stereo descriptors
    # In a steroid skeleton, the carbon 3 is typically in the A-ring, next to a potentially aromatic ring.
    hydroxy_beta_pattern = Chem.MolFromSmarts('[C@H](O)[C@@H]1CC2')
    if not mol.HasSubstructMatch(hydroxy_beta_pattern):
        return False, "3beta-hydroxy group not properly oriented or absent"

    return True, "Molecule correctly classified as a 3beta-hydroxy steroid"