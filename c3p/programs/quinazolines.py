"""
Classifies: CHEBI:38530 quinazolines
"""
from rdkit import Chem

def is_quinazolines(smiles: str):
    """
    Determines if a molecule is a quinazoline or a substituted derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinazoline or substituted derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the quinazoline core structure
    quinazoline_smarts = "c1ncc2c(n1)cccc2"
    quinazoline_core = Chem.MolFromSmarts(quinazoline_smarts)

    if quinazoline_core is None:
        return False, "Invalid SMARTS pattern for quinazoline core"

    # Check if the molecule contains the quinazoline core
    if mol.HasSubstructMatch(quinazoline_core):
        return True, "Contains quinazoline core"
    else:
        return False, "Does not contain quinazoline core"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38530',
                          'name': 'quinazolines',
                          'definition': 'Any organic heterobicyclic compound '
                                        'based on a quinazoline skeleton and '
                                        'its substituted derivatives.',
                          'parents': ['CHEBI:27171', 'CHEBI:38101']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 50,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 17,
    'precision': 1.0,
    'recall': 0.746268656716418,
    'f1': 0.8547008547008548,
    'accuracy': None}