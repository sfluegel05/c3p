"""
Classifies: CHEBI:88061 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine (contains two or more amino groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all amino groups (NH2, NH, N)
    amino_groups = []
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'N':
            num_hydrogens = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'H')
            if num_hydrogens > 0 or (atom.GetFormalCharge() == 0 and atom.GetDegree() > 1):
                amino_groups.append(atom)

    if len(amino_groups) >= 2:
        return True, f"Contains {len(amino_groups)} amino groups"
    else:
        return False, f"Contains only {len(amino_groups)} amino group(s)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:88061',
                          'name': 'polyamine',
                          'definition': 'Any organic amino compound that '
                                        'contains two or more amino groups.',
                          'parents': ['CHEBI:50047']},
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
    'num_true_positives': 8,
    'num_false_positives': 10,
    'num_true_negatives': 7,
    'num_false_negatives': 9,
    'precision': 0.4444444444444444,
    'recall': 0.47058823529411764,
    'f1': 0.45714285714285713,
    'accuracy': None}