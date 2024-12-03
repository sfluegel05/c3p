"""
Classifies: CHEBI:22798 beta-D-glucoside
"""
from rdkit import Chem

def is_beta_D_glucoside(smiles: str):
    """
    Determines if a molecule is a beta-D-glucoside.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-glucoside substructure with beta-configuration at the anomeric carbon
    beta_D_glucoside = Chem.MolFromSmarts('OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O')

    if beta_D_glucoside is None:
        return False, "Unable to define beta-D-glucoside substructure"

    if mol.HasSubstructMatch(beta_D_glucoside):
        return True, "Contains beta-D-glucoside substructure"
    else:
        return False, "Does not contain beta-D-glucoside substructure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22798',
                          'name': 'beta-D-glucoside',
                          'definition': 'Any D-glucoside in which the anomeric '
                                        'centre has beta-configuration.',
                          'parents': ['CHEBI:35436', 'CHEBI:60980']},
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
    'num_true_positives': 74,
    'num_false_positives': 2,
    'num_true_negatives': 18,
    'num_false_negatives': 0,
    'precision': 0.9736842105263158,
    'recall': 1.0,
    'f1': 0.9866666666666666,
    'accuracy': None}