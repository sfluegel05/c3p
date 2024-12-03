"""
Classifies: CHEBI:25750 oxime
"""
from rdkit import Chem

def is_oxime(smiles: str):
    """
    Determines if a molecule is an oxime (R2C=NOH derived from aldehydes or ketones with hydroxylamine).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxime, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for oxime functionality
    oxime_pattern = Chem.MolFromSmarts("[C]=[N][OH]")

    if mol.HasSubstructMatch(oxime_pattern):
        return True, "Molecule contains oxime functionality"

    return False, "Molecule does not contain oxime functionality"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25750',
                          'name': 'oxime',
                          'definition': 'Compounds of structure R2C=NOH '
                                        'derived from condensation of '
                                        'aldehydes or ketones with '
                                        'hydroxylamine. Oximes from aldehydes '
                                        'may be called aldoximes; those from '
                                        'ketones may be called ketoximes.',
                          'parents': ['CHEBI:50860', 'CHEBI:51143']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 9,
    'num_false_positives': 0,
    'num_true_negatives': 11,
    'num_false_negatives': 2,
    'precision': 1.0,
    'recall': 0.8181818181818182,
    'f1': 0.9,
    'accuracy': None}