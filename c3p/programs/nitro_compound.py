"""
Classifies: CHEBI:35715 nitro compound
"""
from rdkit import Chem

def is_nitro_compound(smiles: str):
    """
    Determines if a molecule is a nitro compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitro compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the nitro group SMARTS pattern
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')
    
    # Check if the molecule contains the nitro group
    if mol.HasSubstructMatch(nitro_pattern):
        return True, "Contains nitro group"
    else:
        return False, "Does not contain nitro group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35715',
                          'name': 'nitro compound',
                          'definition': 'A compound having a nitro group, -NO2 '
                                        '(free valence on nitrogen), which may '
                                        'be attached to carbon, nitrogen (as '
                                        'in nitramines), or oxygen (as in '
                                        'nitrates), among other elements (in '
                                        'the absence of specification, C-nitro '
                                        'compounds are usually implied).',
                          'parents': ['CHEBI:51143']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 22-23: malformed \\N character escape (<string>, line '
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