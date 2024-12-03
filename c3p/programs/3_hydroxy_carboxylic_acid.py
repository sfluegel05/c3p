"""
Classifies: CHEBI:61355 3-hydroxy carboxylic acid
"""
from rdkit import Chem

def is_3_hydroxy_carboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-hydroxy carboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    carboxylic_acid = None
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            neighbors = atom.GetNeighbors()
            oxygens = [n for n in neighbors if n.GetAtomicNum() == 8]
            if len(oxygens) == 2:
                if any(o.GetTotalNumHs() == 1 for o in oxygens):
                    carboxylic_acid = atom
                    break
    if carboxylic_acid is None:
        return False, "No carboxylic acid group found"

    carboxylic_acid_neighbors = carboxylic_acid.GetNeighbors()
    for neighbor in carboxylic_acid_neighbors:
        if neighbor.GetAtomicNum() == 6:  # Carbon
            alpha_carbon = neighbor
            alpha_neighbors = alpha_carbon.GetNeighbors()
            for alpha_neighbor in alpha_neighbors:
                if alpha_neighbor.GetAtomicNum() == 6 and alpha_neighbor != carboxylic_acid:
                    beta_carbon = alpha_neighbor
                    beta_neighbors = beta_carbon.GetNeighbors()
                    for beta_neighbor in beta_neighbors:
                        if beta_neighbor.GetAtomicNum() == 8:  # Oxygen
                            if beta_neighbor.GetTotalNumHs() == 1:
                                return True, "3-hydroxy carboxylic acid found"
    return False, "No hydroxy group beta- to the carboxylic acid group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:61355',
                          'name': '3-hydroxy carboxylic acid',
                          'definition': 'Any hydroxy carboxylic acid which '
                                        'contains a hydroxy group located '
                                        'beta- to the carboxylic acid group.',
                          'parents': ['CHEBI:24669']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 22,
    'num_false_positives': 20,
    'num_true_negatives': 0,
    'num_false_negatives': 4,
    'precision': 0.5238095238095238,
    'recall': 0.8461538461538461,
    'f1': 0.6470588235294118,
    'accuracy': None}