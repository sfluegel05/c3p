"""
Classifies: CHEBI:24397 glycophospholipid
"""
from rdkit import Chem

def is_glycophospholipid(smiles: str):
    """
    Determines if a molecule is a glycophospholipid, defined as any phospholipid that contain both phosphate and carbohydrate as integral structural components.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycophospholipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphate group
    phosphate = Chem.MolFromSmarts('P(=O)(O)(O)O')
    if not mol.HasSubstructMatch(phosphate):
        return False, "No phosphate group found"

    # Check for carbohydrate group
    carbohydrate_patterns = [
        Chem.MolFromSmarts('OC[C@H](O)[C@@H](O)[C@H](O)[C@H](O)C=O'), # common carbohydrate pattern
        Chem.MolFromSmarts('OC[C@H](O)[C@H](O)[C@H](O)[C@@H](O)CO'), # another common pattern
        Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@H](O)[C@H](O)[C@@H](O)[C@H]1O'), # cyclic carbohydrate
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in carbohydrate_patterns):
        return False, "No carbohydrate group found"

    return True, "Molecule contains both phosphate and carbohydrate groups"

# Example usage:
# smiles = "P(O[C@H]1[C@H](O[C@H]2OC([C@@H](O)C(O)[C@H]2O)CO)C(O)[C@H](O)C(O)C1O)(OC[C@H](OC(=O)CCCCCCC/C=C\C/C=C\CCC)COC(=O)CCCCCCC/C=C\CCCCCCCC)(O)=O"
# result, reason = is_glycophospholipid(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24397',
                          'name': 'glycophospholipid',
                          'definition': 'Any phospholipid that contain both '
                                        'phosphate and carbohydrate as '
                                        'integral structural components.',
                          'parents': ['CHEBI:16247', 'CHEBI:33563']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 16,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 5,
    'precision': 1.0,
    'recall': 0.7619047619047619,
    'f1': 0.8648648648648648,
    'accuracy': None}