"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylate group [C(=O)[O-]]
    carboxylate = Chem.MolFromSmarts('C(=O)[O-]')
    carboxylate_matches = mol.GetSubstructMatches(carboxylate)
    if not carboxylate_matches:
        return False, "No carboxylate group found"

    # Check for 2-oxo group [C(=O)C]
    oxo_group = Chem.MolFromSmarts('C(=O)C')
    oxo_matches = mol.GetSubstructMatches(oxo_group)
    if not oxo_matches:
        return False, "No 2-oxo group found"

    # Verify that the oxo group is at the 2-position relative to the carboxylate
    for carboxylate_match in carboxylate_matches:
        carboxylate_carbon = carboxylate_match[0]
        for oxo_match in oxo_matches:
            if carboxylate_carbon in oxo_match and (oxo_match[0] == carboxylate_carbon or oxo_match[1] == carboxylate_carbon):
                return True, "2-oxo monocarboxylic acid anion found"

    return False, "2-oxo group is not at the 2-position relative to the carboxylate"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35179',
                          'name': '2-oxo monocarboxylic acid anion',
                          'definition': 'An oxo monocarboxylic acid anion in '
                                        'which the oxo group is located at the '
                                        '2-position.',
                          'parents': ['CHEBI:35902']},
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
    'num_true_positives': 10,
    'num_false_positives': 0,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}