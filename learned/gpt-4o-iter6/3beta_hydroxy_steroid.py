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

    # Look for steroid tetracyclic backbone (3 six-membered rings followed by one five-membered ring)
    steroid_pattern = Chem.MolFromSmarts('C1C2CC3CC4CCCC(C4)C3C(C2)C1')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"

    # Look for 3beta-hydroxy group (OH group bound with specific stereochemistry)
    # Note that 'O[C@H]' or 'O[C@@H]' depends on the stereochemistry for the beta orientation
    hydroxy_pattern = Chem.MolFromSmarts('C[C@H](O)C')
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "3beta-hydroxy group not found or not in beta orientation"

    return True, "Molecule classified as a 3beta-hydroxy steroid"

# This function takes a SMILES string and checks for the presence of a steroid backbone
# along with a 3beta-oriented hydroxyl group attached to it. Adjustments to stereochemistry 
# in hydroxy pattern might be needed based on actual chemical definition of beta orientation.