"""
Classifies: CHEBI:35753 tricarboxylic acid anion
"""
from rdkit import Chem

def is_tricarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a tricarboxylic acid anion.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tricarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    carboxylate_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = [neighbor.GetSymbol() for neighbor in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and any(neighbor.GetFormalCharge() == -1 for neighbor in atom.GetNeighbors()):
                carboxylate_count += 1

    if carboxylate_count >= 3:
        return True, "Molecule contains at least three carboxylate groups"
    
    return False, "Molecule does not contain at least three carboxylate groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35753',
                          'name': 'tricarboxylic acid anion',
                          'definition': 'Any anion of a tricarboxylic acid  '
                                        'formed by deprotonation of at least '
                                        'one carboxy group.',
                          'parents': ['CHEBI:29067']},
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
    'num_true_positives': 12,
    'num_false_positives': 0,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 1.0,
    'recall': 1.0,
    'f1': 1.0,
    'accuracy': None}