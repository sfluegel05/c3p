"""
Classifies: CHEBI:2468 secondary alpha-hydroxy ketone
"""
from rdkit import Chem

def is_secondary_alpha_hydroxy_ketone(smiles: str):
    """
    Determines if a molecule is a secondary alpha-hydroxy ketone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary alpha-hydroxy ketone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of ketone group
    ketone_pattern = Chem.MolFromSmarts('C(=O)')
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"

    # Check for presence of alpha-hydroxy group
    alpha_hydroxy_pattern = Chem.MolFromSmarts('C(O)[C](=O)')
    if not mol.HasSubstructMatch(alpha_hydroxy_pattern):
        return False, "No alpha-hydroxy group found"

    # Check for secondary alpha-hydroxy ketone
    secondary_alpha_hydroxy_ketone_pattern = Chem.MolFromSmarts('C(O)[C](=O)C')
    if not mol.HasSubstructMatch(secondary_alpha_hydroxy_ketone_pattern):
        return False, "Not a secondary alpha-hydroxy ketone"

    # Ensure the alpha carbon has one hydrogen and one organyl group
    for match in mol.GetSubstructMatches(secondary_alpha_hydroxy_ketone_pattern):
        alpha_carbon = mol.GetAtomWithIdx(match[0])
        if alpha_carbon.GetDegree() == 3 and alpha_carbon.GetTotalNumHs() == 1:
            return True, "Molecule is a secondary alpha-hydroxy ketone"

    return False, "Alpha carbon does not have the required substituents"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:2468',
                          'name': 'secondary alpha-hydroxy ketone',
                          'definition': 'An alpha-hydroxy ketone in which the '
                                        'carbonyl group and the hydroxy group '
                                        'are linked by a carbon bearing one '
                                        'hydrogen and one organyl group. '
                                        'Secondary alpha-hydroxy ketones are '
                                        'also known as acyloins, and are '
                                        'formally derived from reductive '
                                        'coupling of two carboxylic acid '
                                        'groups.',
                          'parents': ['CHEBI:139588', 'CHEBI:35681']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 14,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 10,
    'precision': 1.0,
    'recall': 0.5833333333333334,
    'f1': 0.7368421052631579,
    'accuracy': None}