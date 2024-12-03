"""
Classifies: CHEBI:21731 N-glycosyl compound
"""
from rdkit import Chem

def is_N_glycosyl_compound(smiles: str):
    """
    Determines if a molecule is an N-glycosyl compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-glycosyl compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the substructure pattern for N-glycosyl bond
    n_glycosyl_pattern = Chem.MolFromSmarts("[N;!R][C;!R]1[O;!R][C;!R][C;!R][C;!R][O;!R]1")

    if mol.HasSubstructMatch(n_glycosyl_pattern):
        return True, "N-glycosyl compound detected"
    
    return False, "No N-glycosyl bond detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:21731',
                          'name': 'N-glycosyl compound',
                          'definition': 'A glycosyl compound arising formally '
                                        'from the elimination of water from a '
                                        'glycosidic hydroxy group and an H '
                                        'atom bound to a nitrogen atom, thus '
                                        'creating a C-N bond.',
                          'parents': [   'CHEBI:35352',
                                         'CHEBI:63161',
                                         'CHEBI:63299']},
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
    'num_false_negatives': 71,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}