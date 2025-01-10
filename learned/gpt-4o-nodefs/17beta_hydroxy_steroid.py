"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid contains a steroid backbone with a hydroxyl group at the 17-beta position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for steroid skeleton (cyclopenta[a]phenanthrene)
    steroid_skeleton_pattern = Chem.MolFromSmarts("C1CC2CCC3C4C=CC(=O)CC4CC3C2(C)C1")
    if not mol.HasSubstructMatch(steroid_skeleton_pattern):
        return False, "No steroid skeleton found"
    
    # Define SMARTS pattern for hydroxyl group at 17-beta position
    # This pattern assumes the steroid skeleton has been identified in its correct orientation
    hydroxyl_17beta_pattern = Chem.MolFromSmarts("[C@H]1...O")  # Simplified representation
    if not mol.HasSubstructMatch(hydroxyl_17beta_pattern):
        return False, "No hydroxyl group at 17-beta position found"
    
    return True, "Molecule matches 17beta-hydroxy steroid structure"

# Added metadata based on the chemical class and recognition mechanism
__metadata__ = {   
    'chemical_class': {   
        'name': '17beta-hydroxy steroid',
        'definition': 'Steroids with 17-beta hydroxyl group',
    },
    'message': None,
    'success': True,
    'f1_threshold': 0.8,
    'max_attempts': 5,
}