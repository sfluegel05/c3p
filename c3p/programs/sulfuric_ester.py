"""
Classifies: CHEBI:26819 sulfuric ester
"""
from rdkit import Chem

def is_sulfuric_ester(smiles: str):
    """
    Determines if a molecule is a sulfuric ester (an ester of an alcohol and sulfuric acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfuric ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the sulfuric acid ester group (OS(=O)(=O)O)
    sulfuric_ester_pattern = Chem.MolFromSmarts('OS(=O)(=O)O')
    if mol.HasSubstructMatch(sulfuric_ester_pattern):
        return True, "Contains sulfuric ester group"
    else:
        return False, "Does not contain sulfuric ester group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26819',
                          'name': 'sulfuric ester',
                          'definition': 'An ester of an alcohol and sulfuric '
                                        'acid.',
                          'parents': ['CHEBI:35701', 'CHEBI:37826']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 124-125: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}