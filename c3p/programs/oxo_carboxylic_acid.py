"""
Classifies: CHEBI:25754 oxo carboxylic acid
"""
from rdkit import Chem

def is_oxo_carboxylic_acid(smiles: str):
    """
    Determines if a molecule is an oxo carboxylic acid (contains both an aldehyde/ketone group and a carboxylic acid group).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an oxo carboxylic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if neighbors.count('O') == 2 and neighbors.count('H') == 1:
                carboxylic_acid = True
                break

    if not carboxylic_acid:
        return False, "No carboxylic acid group found"

    # Check for aldehyde (CHO) or ketone (C=O) group
    aldehyde_or_ketone = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if neighbors.count('O') == 1 and neighbors.count('H') == 1:
                aldehyde_or_ketone = True
                break
        elif atom.GetSymbol() == 'C' and atom.GetTotalDegree() == 3:
            neighbors = [nbr.GetSymbol() for nbr in atom.GetNeighbors()]
            if neighbors.count('O') == 1 and neighbors.count('C') == 2:
                aldehyde_or_ketone = True
                break

    if not aldehyde_or_ketone:
        return False, "No aldehyde or ketone group found"

    return True, "Molecule contains both a carboxylic acid group and an aldehyde/ketone group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25754',
                          'name': 'oxo carboxylic acid',
                          'definition': 'Any compound that has an aldehydic or '
                                        'ketonic group as well as a carboxylic '
                                        'acid group in the same molecule.',
                          'parents': ['CHEBI:33575']},
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