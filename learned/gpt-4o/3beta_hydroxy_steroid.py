"""
Classifies: CHEBI:36836 3beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 3beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define flexible steroid backbone pattern
    # Allow for various stereochemistry at each junction - using '@' for stereocenters
    steroid_pattern = Chem.MolFromSmarts('[#6]12[#6][#6][#6]3[#6]1[#6][#6][#6]4[#6]3[#6][#6][#6][#6]2')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone found"
        
    # Define 3beta-hydroxy pattern - using '@' for flexibility in stereocenter recognition
    hydroxy_3beta_pattern = Chem.MolFromSmarts('[C@H]1(O)[#6][#6][#6]2[#6]1[#6][#6][#6]3[#6]2[#6][#6][#6]4[#6]3[#6][#6][#6]1')
    if not mol.HasSubstructMatch(hydroxy_3beta_pattern):
        return False, "3beta-hydroxy group not found"
    
    return True, "Contains steroid backbone and a 3beta-hydroxy group"