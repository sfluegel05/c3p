"""
Classifies: CHEBI:24669 hydroxy carboxylic acid
"""
from rdkit import Chem

def is_hydroxy_carboxylic_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy carboxylic acid (any carboxylic acid with at least one hydroxy group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy carboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    carboxylic_acid = False
    hydroxy_group = False

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            # Check for carboxylic acid group (C(=O)O)
            if neighbors.count('O') == 2 and any(neighbor.GetSymbol() == 'O' and neighbor.GetFormalCharge() == 0 for neighbor in atom.GetNeighbors()):
                carboxylic_acid = True
        if atom.GetSymbol() == 'O':
            # Check for hydroxy group (O-H)
            if any(neighbor.GetSymbol() == 'H' for neighbor in atom.GetNeighbors()):
                hydroxy_group = True
            # Check for hydroxy group (O-C)
            if any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
                hydroxy_group = True

    if carboxylic_acid and hydroxy_group:
        return True, "Molecule is a hydroxy carboxylic acid"
    elif not carboxylic_acid:
        return False, "No carboxylic acid group found"
    elif not hydroxy_group:
        return False, "No hydroxy group found"

    return None, None


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24669',
                          'name': 'hydroxy carboxylic acid',
                          'definition': 'Any carboxylic acid with at least one '
                                        'hydroxy group.',
                          'parents': ['CHEBI:33575', 'CHEBI:33822']},
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
    'num_true_positives': 130,
    'num_false_positives': 10,
    'num_true_negatives': 10,
    'num_false_negatives': 0,
    'precision': 0.9285714285714286,
    'recall': 1.0,
    'f1': 0.962962962962963,
    'accuracy': None}