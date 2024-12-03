"""
Classifies: CHEBI:51614 diarylmethane
"""
from rdkit import Chem

def is_diarylmethane(smiles: str):
    """
    Determines if a molecule is a diarylmethane (two aryl groups connected by a single carbon atom).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diarylmethane, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a carbon atom connected to exactly two aromatic rings
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = atom.GetNeighbors()
            aromatic_neighbors = [neighbor for neighbor in neighbors if neighbor.GetIsAromatic()]
            if len(aromatic_neighbors) == 2:
                return True, "Contains two aryl groups connected by a single carbon atom"

    return False, "Does not contain two aryl groups connected by a single carbon atom"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51614',
                          'name': 'diarylmethane',
                          'definition': 'Any compound containing two aryl '
                                        'groups connected by a single C atom.',
                          'parents': ['CHEBI:33659']},
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
    'num_true_positives': 37,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 0,
    'precision': 0.7708333333333334,
    'recall': 1.0,
    'f1': 0.8705882352941177,
    'accuracy': None}