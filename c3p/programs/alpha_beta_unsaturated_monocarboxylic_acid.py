"""
Classifies: CHEBI:79020 alpha,beta-unsaturated monocarboxylic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_beta_unsaturated_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is an alpha,beta-unsaturated monocarboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha,beta-unsaturated monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    if not mol.HasSubstructMatch(carboxylic_acid):
        return False, "No carboxylic acid group found"

    # Check for alpha,beta-unsaturation (C=C or C#C adjacent to carboxylic acid)
    alpha_beta_unsaturation = Chem.MolFromSmarts('C(=O)C=C')
    if mol.HasSubstructMatch(alpha_beta_unsaturation):
        return True, "Alpha,beta-unsaturated monocarboxylic acid with C=C bond found"

    alpha_beta_unsaturation = Chem.MolFromSmarts('C(=O)C#C')
    if mol.HasSubstructMatch(alpha_beta_unsaturation):
        return True, "Alpha,beta-unsaturated monocarboxylic acid with C#C bond found"

    return False, "No alpha,beta-unsaturation found adjacent to carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:79020',
                          'name': 'alpha,beta-unsaturated monocarboxylic acid',
                          'definition': 'A monocarboxylic acid in which the '
                                        'carbon of the carboxy group is '
                                        'directly attached to a C=C or C#C '
                                        'bond.',
                          'parents': ['CHEBI:25384']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 56,
    'num_false_positives': 3,
    'num_true_negatives': 17,
    'num_false_negatives': 8,
    'precision': 0.9491525423728814,
    'recall': 0.875,
    'f1': 0.9105691056910569,
    'accuracy': None}