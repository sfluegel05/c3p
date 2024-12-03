"""
Classifies: CHEBI:35692 dicarboxylic acid
"""
from rdkit import Chem

def is_dicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a dicarboxylic acid (contains two carboxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a dicarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    carboxy_groups = 0

    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 3:
                oxygen_count = sum(1 for neighbor in neighbors if neighbor.GetSymbol() == 'O')
                hydrogen_count = sum(1 for neighbor in neighbors if neighbor.GetSymbol() == 'H')
                if oxygen_count == 2 and hydrogen_count == 0:
                    carboxy_groups += 1

    if carboxy_groups == 2:
        return True, "Molecule contains two carboxy groups"
    else:
        return False, f"Molecule contains {carboxy_groups} carboxy groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35692',
                          'name': 'dicarboxylic acid',
                          'definition': 'Any carboxylic acid containing two '
                                        'carboxy groups.',
                          'parents': ['CHEBI:131927', 'CHEBI:33575']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 10-11: malformed \\N character escape (<string>, line '
             '1)',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}