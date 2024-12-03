"""
Classifies: CHEBI:24586 hexosamine
"""
from rdkit import Chem

def is_hexosamine(smiles: str):
    """
    Determines if a molecule is a hexosamine (Any 6-carbon amino monosaccharide with at least one alcoholic hydroxy group replaced by an amino group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexosamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for 6 carbon atoms in the main chain
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count < 6:
        return False, "Molecule does not have 6 carbon atoms"

    # Check for at least one amino group
    amino_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            amino_group = True
            break

    if not amino_group:
        return False, "No amino group found"

    # Check for at least one hydroxy group
    hydroxy_group = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and any(neighbor.GetSymbol() == 'C' for neighbor in atom.GetNeighbors()):
            hydroxy_group = True
            break

    if not hydroxy_group:
        return False, "No hydroxy group found"

    # Check if at least one hydroxy group is replaced by an amino group
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            if 'O' in neighbors and 'N' in neighbors:
                return True, "Hexosamine with an amino group replacing a hydroxy group found"

    return False, "No hydroxy group replaced by an amino group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24586',
                          'name': 'hexosamine',
                          'definition': 'Any 6-carbon amino monosaccharide '
                                        'with at least one alcoholic hydroxy '
                                        'group replaced by an amino group.',
                          'parents': ['CHEBI:60926']},
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
    'num_true_positives': 20,
    'num_false_positives': 11,
    'num_true_negatives': 9,
    'num_false_negatives': 5,
    'precision': 0.6451612903225806,
    'recall': 0.8,
    'f1': 0.7142857142857142,
    'accuracy': None}