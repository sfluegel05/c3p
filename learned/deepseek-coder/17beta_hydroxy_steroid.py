"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: CHEBI:28800 17beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17beta-hydroxy steroid based on its SMILES string.
    A 17beta-hydroxy steroid has a steroid backbone with a hydroxyl group at the 17th position in the beta-configuration.

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

    # Define the steroid backbone pattern (four fused rings)
    steroid_backbone_pattern = Chem.MolFromSmarts("[C@]12[C@]3([C@]4([C@]5([C@]1([C@H]([C@H]2[C@H]3[C@H]4[C@H]5)[H])[H])[H])[H])[H]")
    if not mol.HasSubstructMatch(steroid_backbone_pattern):
        return False, "No steroid backbone found"

    # Define the pattern for a hydroxyl group at the 17th position in beta-configuration
    # The beta-configuration means the hydroxyl group is on the same side as the angular methyl group at position 18
    beta_hydroxy_17_pattern = Chem.MolFromSmarts("[C@]12[C@]3([C@]4([C@]5([C@]1([C@H]([C@H]2[C@H]3[C@H]4[C@H]5)[H])[H])[H])[H])[C@@H](O)[C@]1([H])CC[C@@H]1C")
    if not mol.HasSubstructMatch(beta_hydroxy_17_pattern):
        return False, "No 17beta-hydroxy group found"

    return True, "Contains a steroid backbone with a 17beta-hydroxy group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28800',
                          'name': '17beta-hydroxy steroid',
                          'definition': 'A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.',
                          'parents': ['CHEBI:28800', 'CHEBI:28800']},
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