"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find carboxylic acid group (COOH)
    carboxyl_group = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # carbon
            neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
            if neighbors.count(8) == 2 and neighbors.count(1) == 1:
                carboxyl_group = atom
                break

    if carboxyl_group is None:
        return False, "No carboxylic acid group found"

    # Find alpha carbon (carbon attached to carboxyl group)
    alpha_carbon = None
    for neighbor in carboxyl_group.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # carbon
            alpha_carbon = neighbor
            break

    if alpha_carbon is None:
        return False, "No alpha carbon found"

    # Find amino group (NH2) attached to alpha carbon
    amino_group = None
    for neighbor in alpha_carbon.GetNeighbors():
        if neighbor.GetAtomicNum() == 7:  # nitrogen
            if sum(1 for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 1) == 2:
                amino_group = neighbor
                break

    if amino_group is None:
        return False, "No amino group found on alpha carbon"

    return True, "Alpha-amino acid detected"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33704',
                          'name': 'alpha-amino acid',
                          'definition': 'An amino acid in which the amino '
                                        'group is located on the carbon atom '
                                        'at the position alpha to the carboxy '
                                        'group.',
                          'parents': ['CHEBI:33709']},
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
    'num_false_negatives': 107,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}