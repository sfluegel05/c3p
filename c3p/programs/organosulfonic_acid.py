"""
Classifies: CHEBI:33551 organosulfonic acid
"""
from rdkit import Chem

def is_organosulfonic_acid(smiles: str):
    """
    Determines if a molecule is an organosulfonic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organosulfonic acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the sulfonic acid group pattern
    sulfonic_acid_pattern = Chem.MolFromSmarts('S(=O)(=O)[O]')
    
    # Check if the molecule contains a sulfonic acid group
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group found"

    # Get all matches for the sulfonic acid group
    matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
    
    for match in matches:
        sulfonic_acid_atom = match[0]
        # Check if the sulfur atom is directly connected to a carbon atom
        sulfur_atom = mol.GetAtomWithIdx(sulfonic_acid_atom)
        neighbors = sulfur_atom.GetNeighbors()
        if any(neighbor.GetSymbol() == 'C' for neighbor in neighbors):
            return True, "Sulfonic acid group is directly linked to a carbon atom"
    
    return False, "Sulfonic acid group is not directly linked to a carbon atom"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33551',
                          'name': 'organosulfonic acid',
                          'definition': 'An organic derivative of sulfonic '
                                        'acid in which the sulfo group is '
                                        'linked directly to carbon.',
                          'parents': [   'CHEBI:33261',
                                         'CHEBI:33552',
                                         'CHEBI:64709']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 54-55: malformed \\N character escape (<string>, line '
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