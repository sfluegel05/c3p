"""
Classifies: CHEBI:35972 dihydroxy monocarboxylic acid
"""
from rdkit import Chem

def is_dihydroxy_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dihydroxy monocarboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dihydroxy monocarboxylic acid, False otherwise
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
            if neighbors.count('O') == 2 and 'C' in neighbors:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for at least two hydroxy groups (OH)
    hydroxy_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(n.GetSymbol() == 'H' for n in atom.GetNeighbors()))

    if hydroxy_count < 2:
        return False, "Less than two hydroxy groups found"

    return True, "Molecule is a dihydroxy monocarboxylic acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35972',
                          'name': 'dihydroxy monocarboxylic acid',
                          'definition': 'Any  hydroxy monocarboxylic acid '
                                        'carrying at least two hydroxy groups.',
                          'parents': ['CHEBI:35868']},
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
    'num_true_negatives': 14,
    'num_false_negatives': 14,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}