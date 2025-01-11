"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
"""
Classifies: CHEBI:28875 myo-inositol phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate with a myo-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the myo-inositol backbone pattern (6-membered ring with hydroxyl or phosphate groups)
    myo_inositol_pattern = Chem.MolFromSmarts("[C@H]1([C@@H]([OH,OP])([C@H]([OH,OP])([C@@H]([OH,OP])([C@H]([OH,OP])([C@@H]1[OH,OP])))))")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol backbone found"

    # Look for at least one phosphate group (-OP(O)(O)=O or -OP(O)(O)-)
    phosphate_pattern = Chem.MolFromSmarts("[OX2][PX4](=[OX1])([OX2H,OX2-])[OX2H,OX2-]")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"

    return True, "Contains myo-inositol backbone with at least one phosphate group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28875',
                          'name': 'myo-inositol phosphate',
                          'definition': 'An inositol phosphate in which the inositol component has myo-configuration.',
                          'parents': ['CHEBI:28874', 'CHEBI:24838']},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                  'f1_threshold': 0.8,
                  'max_attempts': 5,
                  'max_positive_instances': None,
                  'max_positive_to_test': None,
                  'max_negative_to_test': None,
                  'max_positive_in_prompt': 50,
                  'max_negative_in_prompt': 20,
                  'max_instances_in_prompt': 100,
                  'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 150,
    'num_false_positives': 4,
    'num_true_negatives': 182407,
    'num_false_negatives': 23,
    'num_negatives': None,
    'precision': 0.974025974025974,
    'recall': 0.8670520231213873,
    'f1': 0.9174311926605504,
    'accuracy': 0.9998521228585199}