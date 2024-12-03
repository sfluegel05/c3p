"""
Classifies: CHEBI:37485 hydroxyanthraquinones
"""
from rdkit import Chem

def is_hydroxyanthraquinones(smiles: str):
    """
    Determines if a molecule is a hydroxyanthraquinone (an anthraquinone substituted with one or more hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxyanthraquinone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define anthraquinone core pattern
    anthraquinone_smarts = 'O=C1C2=CC=CC=C2C(=O)C3=CC=CC=C13'
    anthraquinone_core = Chem.MolFromSmarts(anthraquinone_smarts)
    
    if not mol.HasSubstructMatch(anthraquinone_core):
        return False, "No anthraquinone core found"

    # Check for hydroxy groups attached to the anthraquinone core
    hydroxy_smarts = '[OH]'
    hydroxy_groups = mol.GetSubstructMatches(Chem.MolFromSmarts(hydroxy_smarts))
    
    if not hydroxy_groups:
        return False, "No hydroxy groups found"

    return True, "Hydroxyanthraquinone found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37485',
                          'name': 'hydroxyanthraquinones',
                          'definition': 'An anthraquinone substituted with one '
                                        'or more hydroxy groups.',
                          'parents': ['CHEBI:22580', 'CHEBI:33853']},
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
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 33,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}