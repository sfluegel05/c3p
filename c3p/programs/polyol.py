"""
Classifies: CHEBI:26191 polyol
"""
from rdkit import Chem

def is_polyol(smiles: str):
    """
    Determines if a molecule is a polyol (a compound that contains two or more hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of hydroxy groups (-OH)
    hydroxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'C':
                    if any(n.GetSymbol() == 'H' for n in neighbor.GetNeighbors()):
                        hydroxy_count += 1
                    elif any(n.GetSymbol() == 'C' for n in neighbor.GetNeighbors()):
                        hydroxy_count += 1

    if hydroxy_count >= 2:
        return True, f"Contains {hydroxy_count} hydroxy groups"
    else:
        return False, f"Contains only {hydroxy_count} hydroxy groups"

# Example usage:
# print(is_polyol("COc1ccc(CC(O)C(=O)c2ccc(OC)cc2O)cc1"))  # Should return (True, "Contains 3 hydroxy groups")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26191',
                          'name': 'polyol',
                          'definition': 'A compound that contains two or more '
                                        'hydroxy groups.',
                          'parents': ['CHEBI:33822']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '[13:39:42] SMILES Parse Error: syntax error while parsing: '
             'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C\n'
             '[13:39:42] SMILES Parse Error: Failed parsing SMILES '
             "'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C' "
             'for input: '
             "'CC[C@@H]1CC[C@@H]2O[C@@]3(O[C@@H](C[C@@H](C)O)[C@@H](C)CC3=O)[C@@H](C)[C@H](OC(=O)\\C=C\\[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@H](C)C(=O)[C@](C)(O)[C@@H](O)[C@H](C)C\\C=C\\C=C\x01)[C@@H]2C'\n",
    'stdout': '',
    'num_true_positives': 147,
    'num_false_positives': 19,
    'num_true_negatives': 1,
    'num_false_negatives': 1,
    'precision': 0.8855421686746988,
    'recall': 0.9932432432432432,
    'f1': 0.9363057324840764,
    'accuracy': None}