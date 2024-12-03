"""
Classifies: CHEBI:35969 3-hydroxy monocarboxylic acid
"""
from rdkit import Chem

def is_3_hydroxy_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy monocarboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups
    carboxylic_acid = Chem.MolFromSmarts('C(=O)O')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid)

    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    for match in carboxylic_acid_matches:
        carboxylic_c = match[0]

        # Check for hydroxy group beta to carboxy group
        for neighbor in mol.GetAtomWithIdx(carboxylic_c).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # carbon
                for beta_neighbor in neighbor.GetNeighbors():
                    if beta_neighbor.GetAtomicNum() == 6:  # carbon
                        for gamma_neighbor in beta_neighbor.GetNeighbors():
                            if gamma_neighbor.GetAtomicNum() == 8:  # oxygen
                                if any([n.GetAtomicNum() == 1 for n in gamma_neighbor.GetNeighbors()]):  # hydrogen
                                    return True, "3-hydroxy monocarboxylic acid found"

    return False, "No 3-hydroxy group beta to carboxy group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35969',
                          'name': '3-hydroxy monocarboxylic acid',
                          'definition': 'A hydroxy monocarboxylic acid that '
                                        'has a hydroxy group beta to the '
                                        'carboxy group.',
                          'parents': ['CHEBI:35868', 'CHEBI:61355']},
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
    'num_true_negatives': 10,
    'num_false_negatives': 10,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}