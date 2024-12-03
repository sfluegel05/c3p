"""
Classifies: CHEBI:35356 dicarboximide
"""
from rdkit import Chem

def is_dicarboximide(smiles: str):
    """
    Determines if a molecule is a dicarboximide (an imide in which the two acyl substituents on nitrogen are carboacyl groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboximide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for dicarboximide
    dicarboximide_smarts = '[NX3](C(=O)*)C(=O)*'
    pattern = Chem.MolFromSmarts(dicarboximide_smarts)

    if mol.HasSubstructMatch(pattern):
        return True, "Molecule matches dicarboximide pattern"
    else:
        return False, "Molecule does not match dicarboximide pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35356',
                          'name': 'dicarboximide',
                          'definition': 'An imide in which the two acyl '
                                        'substituents on nitrogen are '
                                        'carboacyl groups.',
                          'parents': ['CHEBI:24782']},
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
    'num_true_positives': 22,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}