"""
Classifies: CHEBI:17334 penicillin
"""
"""
Classifies: CHEBI:17334 penicillin
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_penicillin(smiles: str):
    """
    Determines if a molecule is a penicillin based on its SMILES string.
    A penicillin is a member of the group of substituted penams containing:
    - Two methyl substituents at position 2
    - A carboxylate substituent at position 3
    - A carboxamido group at position 6

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a penicillin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the penicillin core structure pattern
    penicillin_core = Chem.MolFromSmarts("[C@@H]1[C@H]2N(C1=O)[C@H](C(S2)(C)C)C(=O)O")
    if not mol.HasSubstructMatch(penicillin_core):
        return False, "No penicillin core structure found"

    # Check for two methyl groups at position 2
    methyl_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2N(C1=O)[C@H](C(S2)(C)C)C(=O)O")
    methyl_matches = mol.GetSubstructMatches(methyl_pattern)
    if len(methyl_matches) == 0:
        return False, "No methyl groups found at position 2"

    # Check for carboxylate group at position 3
    carboxylate_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2N(C1=O)[C@H](C(S2)(C)C)C(=O)O")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    if len(carboxylate_matches) == 0:
        return False, "No carboxylate group found at position 3"

    # Check for carboxamido group at position 6
    carboxamido_pattern = Chem.MolFromSmarts("[C@@H]1[C@H]2N(C1=O)[C@H](C(S2)(C)C)C(=O)O")
    carboxamido_matches = mol.GetSubstructMatches(carboxamido_pattern)
    if len(carboxamido_matches) == 0:
        return False, "No carboxamido group found at position 6"

    return True, "Contains penicillin core structure with required substituents"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17334',
                          'name': 'penicillin',
                          'definition': 'Any member of the group of substituted penams containing two methyl substituents at position 2, a carboxylate substituent at position 3 and a carboxamido group at position 6.',
                          'parents': ['CHEBI:17334', 'CHEBI:17334']},
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