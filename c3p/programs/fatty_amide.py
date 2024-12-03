"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide (a monocarboxylic acid amide derived from a fatty acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of the amide group (C(=O)N)
    amide_group = Chem.MolFromSmarts("C(=O)N")
    if not mol.HasSubstructMatch(amide_group):
        return False, "No amide group found"

    # Check for the presence of a long carbon chain characteristic of fatty acids
    carbon_chain = Chem.MolFromSmarts("CCCCCCCC")
    if not mol.HasSubstructMatch(carbon_chain):
        return False, "No long carbon chain found"

    return True, "Fatty amide structure found"

# Example usage:
# print(is_fatty_amide("CCCCCCCCCCCCCCC(=O)NCCO"))  # Should return (True, "Fatty amide structure found")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:29348',
                          'name': 'fatty amide',
                          'definition': 'A monocarboxylic acid amide derived '
                                        'from a fatty acid.',
                          'parents': ['CHEBI:29347', 'CHEBI:61697']},
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
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 17,
    'precision': 1.0,
    'recall': 0.5641025641025641,
    'f1': 0.7213114754098361,
    'accuracy': None}