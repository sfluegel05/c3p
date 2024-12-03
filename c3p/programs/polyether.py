"""
Classifies: CHEBI:46774 polyether
"""
from rdkit import Chem

def is_polyether(smiles: str):
    """
    Determines if a molecule is a polyether (contains more than one ether linkage).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyether, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of ether linkages (C-O-C)
    ether_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            if neighbors.count('C') >= 2:
                ether_count += 1

    if ether_count > 1:
        return True, f"Polyether with {ether_count} ether linkages"
    else:
        return False, "Not a polyether"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46774',
                          'name': 'polyether',
                          'definition': 'Any ether that contains more than one '
                                        'ether linkage.',
                          'parents': ['CHEBI:25698']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 15,
    'num_false_positives': 10,
    'num_true_negatives': 5,
    'num_false_negatives': 0,
    'precision': 0.6,
    'recall': 1.0,
    'f1': 0.7499999999999999,
    'accuracy': None}