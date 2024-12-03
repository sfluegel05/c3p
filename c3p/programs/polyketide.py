"""
Classifies: CHEBI:26188 polyketide
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_polyketide(smiles: str):
    """
    Determines if a molecule is a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyketide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alternating carbonyl and methylene groups
    pattern = Chem.MolFromSmarts('C(=O)-C-C(=O)-C')
    if mol.HasSubstructMatch(pattern):
        return True, "Contains alternating carbonyl and methylene groups"

    return False, "Does not contain alternating carbonyl and methylene groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26188',
                          'name': 'polyketide',
                          'definition': 'Natural and synthetic compounds '
                                        'containing alternating carbonyl and '
                                        'methylene groups '
                                        "('beta-polyketones'), biogenetically "
                                        'derived from repeated condensation of '
                                        'acetyl coenzyme A (via malonyl '
                                        'coenzyme A), and usually the '
                                        'compounds derived from them by '
                                        'further condensations, etc. '
                                        'Considered by many to be synonymous '
                                        'with the less frequently used terms '
                                        'acetogenins and ketides.',
                          'parents': ['CHEBI:36963']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 132-133: malformed \\N character escape (<string>, line '
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