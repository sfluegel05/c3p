"""
Classifies: CHEBI:29347 monocarboxylic acid amide
"""
from rdkit import Chem

def is_monocarboxylic_acid_amide(smiles: str):
    """
    Determines if a molecule is a monocarboxylic acid amide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monocarboxylic acid amide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxamide group (C(=O)N)
    carboxamide_group = Chem.MolFromSmarts('C(=O)N')
    if not mol.HasSubstructMatch(carboxamide_group):
        return False, "No carboxamide group found"

    # Check if the carboxamide group is derived from a monocarboxylic acid
    carboxylic_acid_group = Chem.MolFromSmarts('C(=O)O')
    if mol.HasSubstructMatch(carboxylic_acid_group):
        return False, "Contains a carboxylic acid group"

    # Check if the carboxamide is derived from a monocarboxylic acid
    carboxamide_matches = mol.GetSubstructMatches(carboxamide_group)
    for match in carboxamide_matches:
        carbon = mol.GetAtomWithIdx(match[0])
        if carbon.GetDegree() == 3 and any(neighbor.GetSymbol() == 'O' for neighbor in carbon.GetNeighbors()):
            return True, "Monocarboxylic acid amide found"
    
    return False, "No monocarboxylic acid amide found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:29347',
                          'name': 'monocarboxylic acid amide',
                          'definition': 'A carboxamide derived from a '
                                        'monocarboxylic acid.',
                          'parents': ['CHEBI:37622']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': False,
    'best': True,
    'error': "(unicode error) 'unicodeescape' codec can't decode bytes in "
             'position 11-12: malformed \\N character escape (<string>, line '
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