"""
Classifies: CHEBI:26820 sulfates
"""
from rdkit import Chem

def is_sulfates(smiles: str):
    """
    Determines if a molecule is a sulfate (salt or ester of sulfuric acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfate, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    sulfate_found = False
    reason = "No sulfate group found"

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 16:  # Sulfur atom
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 4:
                oxy_count = sum(1 for neighbor in neighbors if neighbor.GetAtomicNum() == 8)
                if oxy_count == 4:
                    sulfate_found = True
                    reason = "Sulfate group found"
                    break

    return sulfate_found, reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26820',
                          'name': 'sulfates',
                          'definition': 'Salts and esters of sulfuric acid',
                          'parents': ['CHEBI:37826']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 123-124: malformed \\N character escape (<string>, line '
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