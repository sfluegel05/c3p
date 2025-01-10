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
    
    # Define a more general SMARTS pattern for the steroid scaffold (4-ring system)
    steroid_skeleton_pattern = Chem.MolFromSmarts("[R2]1[R2][R3]2[R3][R2][R2]3[R2][R3][R3][R2][R2]4[R3][R2][R2][R2][R2]4[C@R]3[R2]2[R3]1")
    if not mol.HasSubstructMatch(steroid_skeleton_pattern):
        return False, "No steroid skeleton found"
    
    # Define a SMARTS pattern for hydroxyl group at 17-beta position
    # Assuming a broader context of steroid structure, flexible pattern for 17-beta hydroxyl
    hydroxyl_17beta_pattern = Chem.MolFromSmarts("[C@H]1CC[C@]2([C@@H]([R2][R1])CC[C@@]2([H])C)[C@H](O)C=1")
    if not mol.HasSubstructMatch(hydroxyl_17beta_pattern):
        return False, "No hydroxyl group at 17-beta position found"
    
    return True, "Molecule matches 17beta-hydroxy steroid structure"

# Added updated metadata to reflect the flexible recognition mechanisms
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