"""
Classifies: CHEBI:63436 carbohydrate acid derivative
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_carbohydrate_acid_derivative(smiles: str):
    """
    Determines if a molecule is a carbohydrate acid derivative.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbohydrate acid derivative, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of carbohydrate structure
    carbohydrate = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O' and len(atom.GetNeighbors()) == 2:
            neighbors = atom.GetNeighbors()
            if neighbors[0].GetSymbol() == 'C' and neighbors[1].GetSymbol() == 'C':
                carbohydrate = True
                break

    if not carbohydrate:
        return False, "No carbohydrate structure found"

    # Check for presence of carboxylic acid or its derivative
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = atom.GetNeighbors()
            if any(neighbor.GetSymbol() == 'O' and neighbor.GetTotalDegree() == 2 for neighbor in neighbors):
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid or its derivative found"

    return True, "Molecule is a carbohydrate acid derivative"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:63436',
                          'name': 'carbohydrate acid derivative',
                          'definition': 'A carbohydrate derivative that is '
                                        'formally obtained from a carbohydrate '
                                        'acid.',
                          'parents': ['CHEBI:63299']},
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
    'num_true_positives': 60,
    'num_false_positives': 8,
    'num_true_negatives': 12,
    'num_false_negatives': 0,
    'precision': 0.8823529411764706,
    'recall': 1.0,
    'f1': 0.9375,
    'accuracy': None}