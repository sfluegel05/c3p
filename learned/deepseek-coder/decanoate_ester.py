"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:75840 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester is a fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the ester group attached to a 10-carbon chain
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H3]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No decanoate ester group found"

    # Check if the ester is attached to a 10-carbon chain
    decanoate_pattern = Chem.MolFromSmarts("[CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H2][CX4H3]")
    if not mol.HasSubstructMatch(decanoate_pattern):
        return False, "No 10-carbon chain found"

    return True, "Contains a decanoate ester group (10-carbon chain with ester linkage)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:75840',
                          'name': 'decanoate ester',
                          'definition': 'A fatty acid ester resulting from the '
                                        'formal condensation of the carboxy '
                                        'group of decanoic acid (capric acid) '
                                        'with the hydroxy group of an alcohol '
                                        'or phenol.',
                          'parents': ['CHEBI:47778', 'CHEBI:76579']},
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