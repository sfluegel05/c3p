"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
"""
Classifies: CHEBI:38698 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Phore import FragmentMatcher

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid is a steroid with a hydroxyl group attached at the 17th position
    in the alpha configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for steroid scaffold
    steroid_pattern = Chem.MolFromSmarts("[C@]123[C@H]([C@@H]4[C@@]([C@@H]([C@H]([C@@H]2C)C(C)(C)C)C)(C)C=C5C[C@]6([H])C[C@@]7([H])C[C@@H](O)[C@]46C(=O)C=C7[C@]35C)O")
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid scaffold found"
    
    # Check for 17-hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[C@H](O)") 
    hydroxyl_match = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_match) == 0:
        return False, "No hydroxyl group found"
    
    # Check alpha stereochemistry at C17
    frag_matcher = FragmentMatcher()
    frag_matcher.addFragNRing("[C@H](O)[C@@]13[C@H]([C@@H]2[C@@]([C@@H]([C@H]([C@@H]1C)C(C)(C)C)C)(C)C=C4C[C@]5([H])C[C@@]6([H])C[C@]7([H])[C@H]([C@@]56C(=O)C=C7)O4)C(=O)C[C@@H]23", ring_nums=[1,2,3])
    has_alpha_config = frag_matcher.countMatches(mol) > 0
    if not has_alpha_config:
        return False, "Hydroxyl group not in alpha configuration at C17"
    
    # Additional checks for steroid-like properties
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() != 4:
        return False, "Does not have 4 rings"
    
    if sum(1 for ring in ring_info.AtomRings() if ring.IsAromatic()) > 0:
        return False, "Contains aromatic rings, steroids should be fully saturated"
    
    return True, "Steroid with hydroxyl group in alpha configuration at C17"


__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:38698',
        'name': '17alpha-hydroxy steroid',
        'definition': 'The alpha-stereoisomer of 17-hydroxy steroid.',
        'parents': ['CHEBI:38697', 'CHEBI:35796']
    },
    'config': {
        'llm_model_name': 'lbl/claude-sonnet',
        'f1_threshold': 0.8,
        'max_attempts': 5,
        'max_positive_instances': None,
        'max_positive_to_test': None,
        'max_negative_to_test': None,
        'max_positive_in_prompt': 50,
        'max_negative_in_prompt': 20,
        'max_instances_in_prompt': 100,
        'test_proportion': 0.1
    },
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 186,
    'num_false_positives': 0,
    'num_true_negatives': 182421,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9947598253275109,
    'f1': 0.9973544973544974,
    'accuracy': 0.9999945453312859
}