"""
Classifies: CHEBI:26776 stilbenoid
"""
from rdkit import Chem

def is_stilbenoid(smiles: str):
    """
    Determines if a molecule is a stilbenoid (characterized by a 1,2-diphenylethylene backbone).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a stilbenoid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the stilbenoid substructure (1,2-diphenylethylene backbone)
    stilbenoid_substructure = Chem.MolFromSmarts("C=Cc1ccccc1c2ccccc2")

    # Check if the molecule contains the stilbenoid substructure
    if mol.HasSubstructMatch(stilbenoid_substructure):
        return True, "Contains 1,2-diphenylethylene backbone characteristic of stilbenoids"
    else:
        return False, "Does not contain 1,2-diphenylethylene backbone"

# Example usage:
# print(is_stilbenoid("C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C=CC3=CC=CC=C3)O[C@@H]1CN(C)C(=O)NC4CCCCC4)[C@H](C)CO"))


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26776',
                          'name': 'stilbenoid',
                          'definition': 'Any olefinic compound characterised '
                                        'by a 1,2-diphenylethylene backbone.',
                          'parents': ['CHEBI:33659', 'CHEBI:78840']},
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
    'num_true_negatives': 20,
    'num_false_negatives': 56,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': None}