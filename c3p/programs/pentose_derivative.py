"""
Classifies: CHEBI:63409 pentose derivative
"""
from rdkit import Chem

def is_pentose_derivative(smiles: str):
    """
    Determines if a molecule is a pentose derivative (a monosaccharide derivative that is formally obtained from a pentose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pentose derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a pentose (5-carbon sugar) backbone
    pentose_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = atom.GetNeighbors()
            oxygens = [nbr for nbr in neighbors if nbr.GetSymbol() == 'O']
            carbons = [nbr for nbr in neighbors if nbr.GetSymbol() == 'C']
            if len(oxygens) >= 1 and len(carbons) >= 2:
                pentose_found = True
                break

    if pentose_found:
        return True, "Molecule is a pentose derivative"
    else:
        return False, "Molecule is not a pentose derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63409',
                          'name': 'pentose derivative',
                          'definition': 'A monosaccharide derivative that is '
                                        'formally obtained from a pentose.',
                          'parents': ['CHEBI:63367']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 4,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 39,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.6610169491525424,
    'recall': 1.0,
    'f1': 0.7959183673469388,
    'accuracy': None}