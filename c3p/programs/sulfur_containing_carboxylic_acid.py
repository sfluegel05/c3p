"""
Classifies: CHEBI:33576 sulfur-containing carboxylic acid
"""
from rdkit import Chem

def is_sulfur_containing_carboxylic_acid(smiles: str):
    """
    Determines if a molecule is a sulfur-containing carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfur-containing carboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = [n.GetSymbol() for n in atom.GetNeighbors()]
            o_neighbors = [n for n in atom.GetNeighbors() if n.GetSymbol() == 'O']
            if len(o_neighbors) == 2 and any(o.GetTotalValence() == 1 for o in o_neighbors):
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for sulfur substituent
    sulfur_substituent = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'S':
            sulfur_substituent = True
            break

    if not sulfur_substituent:
        return False, "No sulfur substituent found"

    return True, "Sulfur-containing carboxylic acid found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33576',
                          'name': 'sulfur-containing carboxylic acid',
                          'definition': 'Any carboxylic acid having a sulfur '
                                        'substituent.',
                          'parents': ['CHEBI:33261', 'CHEBI:33575']},
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
    'num_true_positives': 1,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 14,
    'precision': 1.0,
    'recall': 0.06666666666666667,
    'f1': 0.125,
    'accuracy': None}