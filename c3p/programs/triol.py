"""
Classifies: CHEBI:27136 triol
"""
from rdkit import Chem

def is_triol(smiles: str):
    """
    Determines if a molecule is a triol (a chemical compound containing three hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of hydroxy groups (OH)
    hydroxy_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            neighbors = atom.GetNeighbors()
            if len(neighbors) == 1 and neighbors[0].GetSymbol() == 'C':
                hydroxy_count += 1

    if hydroxy_count == 3:
        return True, "Molecule contains exactly three hydroxy groups"
    else:
        return False, f"Molecule contains {hydroxy_count} hydroxy groups, not three"

# Example usage:
# smiles = "CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)[C@@](C)(CO)[C@@H]5[C@H](O)C[C@@]34C)[C@@H]2C1)C(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O"
# result, reason = is_triol(smiles)
# print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27136',
                          'name': 'triol',
                          'definition': 'A chemical compound containing three '
                                        'hydroxy groups.',
                          'parents': ['CHEBI:26191']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '',
    'stdout': '',
    'num_true_positives': 9,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 28,
    'precision': 0.9,
    'recall': 0.24324324324324326,
    'f1': 0.3829787234042553,
    'accuracy': None}