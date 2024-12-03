"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a phosphate group
    phosphate_group = Chem.MolFromSmarts('P(=O)(O)(O)O')
    if not mol.HasSubstructMatch(phosphate_group):
        return False, "No phosphate group found"

    # Check for the presence of a monosaccharide
    monosaccharide_group = Chem.MolFromSmarts('[C@H]([O])[C@H]([O])[C@H]([O])[C@H]([O])[C@H]([O])[C@H2]O')
    if not mol.HasSubstructMatch(monosaccharide_group):
        return False, "No monosaccharide structure found"

    # Check for esterification with phosphoric acid
    ester_group = Chem.MolFromSmarts('C-O-P(=O)(O)O')
    if not mol.HasSubstructMatch(ester_group):
        return False, "No esterification with phosphoric acid found"

    return True, "Molecule is a phospho sugar"

# Example usage:
print(is_phospho_sugar('Nc1cncn1C1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O'))  # Example SMILES string


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33447',
                          'name': 'phospho sugar',
                          'definition': 'Any monosaccharide containing an '
                                        'alcoholic hydroxy group esterified '
                                        'with phosphoric acid.',
                          'parents': ['CHEBI:26816', 'CHEBI:63367']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': "(False, 'No monosaccharide structure found')\n",
    'num_true_positives': 5,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 55,
    'precision': 1.0,
    'recall': 0.08333333333333333,
    'f1': 0.15384615384615385,
    'accuracy': None}