"""
Classifies: CHEBI:25529 pyridinecarboxamide
"""
from rdkit import Chem

def is_pyridinecarboxamide(smiles: str):
    """
    Determines if a molecule is a pyridinecarboxamide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyridinecarboxamide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for pyridine ring
    pyridine = Chem.MolFromSmarts('c1ccncc1')
    if not mol.HasSubstructMatch(pyridine):
        return False, "No pyridine ring found"

    # Check for carboxamide group
    carboxamide = Chem.MolFromSmarts('NC(=O)c')
    if not mol.HasSubstructMatch(carboxamide):
        return False, "No carboxamide group found"

    return True, "Molecule is a pyridinecarboxamide"

# Example usage
smiles = "NC(=O)c1ccccn1"
result, reason = is_pyridinecarboxamide(smiles)
print(result, reason)  # Output: True Molecule is a pyridinecarboxamide


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25529',
                          'name': 'pyridinecarboxamide',
                          'definition': 'A member of the class of pyridines '
                                        'that is a substituted pyridine in '
                                        'which at least one of the '
                                        'substituents is a carboxamide or '
                                        'N-substituted caraboxamide group.',
                          'parents': ['CHEBI:26421', 'CHEBI:37622']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'True Molecule is a pyridinecarboxamide\n',
    'num_true_positives': 32,
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 0.8,
    'recall': 1.0,
    'f1': 0.888888888888889,
    'accuracy': None}