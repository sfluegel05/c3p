"""
Classifies: CHEBI:26878 tertiary alcohol
"""
from rdkit import Chem

def is_tertiary_alcohol(smiles: str):
    """
    Determines if a molecule is a tertiary alcohol.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tertiary alcohol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate through all atoms in the molecule
    for atom in mol.GetAtoms():
        # Check if the atom is an oxygen and is part of a hydroxyl group
        if atom.GetSymbol() == 'O' and atom.GetDegree() == 1:
            neighbor = atom.GetNeighbors()[0]
            if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 4:
                # Check if the carbon has exactly three other carbon atoms attached
                carbon_neighbors = [n for n in neighbor.GetNeighbors() if n.GetSymbol() == 'C']
                if len(carbon_neighbors) == 3:
                    return True, "Tertiary alcohol identified"
    
    return False, "No tertiary alcohol structure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26878',
                          'name': 'tertiary alcohol',
                          'definition': 'A tertiary alcohol is a compound in '
                                        'which a hydroxy group, -OH, is '
                                        'attached to a saturated carbon atom '
                                        'which has three other carbon atoms '
                                        'attached to it.',
                          'parents': ['CHEBI:30879']},
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
    'num_true_positives': 106,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 3,
    'precision': 1.0,
    'recall': 0.9724770642201835,
    'f1': 0.986046511627907,
    'accuracy': None}