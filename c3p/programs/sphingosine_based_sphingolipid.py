"""
Classifies: CHEBI:82830 sphingosine-based sphingolipid
"""
from rdkit import Chem


def is_sphingosine_based_sphingolipid(smiles: str):
    """
    Determines if a molecule is a sphingosine-based sphingolipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingosine-based sphingolipid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sphingosine backbone features
    sphingosine_substructure = Chem.MolFromSmarts('C[C@@H](O)[C@H](CO)NC=O')

    if not mol.HasSubstructMatch(sphingosine_substructure):
        return False, "Molecule does not contain the sphingosine backbone"

    return True, "Molecule contains the sphingosine backbone"

# Example usage
smiles_list = [
    "CCCCCCCCCCCCC\C=C\[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCCCC",
    "CCCCCCCCCCCCCCCCCCCCCC[C@@H](O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\C=C\CCCCCCCCCCCCC",
    "CCCCCCCCCCCCC\C=C\[C@@H](O)[C@H](CO)NC(=O)CCCCCCC\C=C/CCCCCCCC",
    "CCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COC1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\C=C\CCCCCCCCCCCCC",
    "CCCCCCCCCCCC\C=C\[C@@H](O)[C@H](COC1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC([*])=O",
    "CCCCCCCCCCCCC\C=C\[C@@H](O)[C@H](CO)NC(C)=O",
    "[C@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCOC(CCCCCCC/C=C\C/C=C\CCCCC)=O)CO",
    "CCCCCCCCCCCC\C=C\[C@@H](O)[C@H](CO)NC([*])=O",
    "CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]1O)[C@H](O)\C=C\CCCCCCCCCCCCC",
    "CCCCCCCCCCCC\C=C\[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]1O)NC(=O)CCCCCCCCCCCCC\C=C/CCCCCCCC",
    "CCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\C=C\CCCCCCCCCCCCC"
]

for smiles in smiles_list:
    result, reason = is_sphingosine_based_sphingolipid(smiles)
    print(f"SMILES: {smiles}, Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:82830',
                          'name': 'sphingosine-based sphingolipid',
                          'definition': 'A family of sphingolipids that share '
                                        'a common structural feature, a '
                                        'sphingosine base backbone.',
                          'parents': ['CHEBI:26739']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCCCC, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n'
              'SMILES: '
              'CCCCCCCCCCCCCCCCCCCCCC[C@@H](O)C(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n'
              'SMILES: '
              'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)CCCCCCC\\C=C/CCCCCCCC, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n'
              'SMILES: '
              'CCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COC1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n'
              'SMILES: '
              'CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](COC1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC([*])=O, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n'
              'SMILES: CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(C)=O, Result: '
              'True, Reason: Molecule contains the sphingosine backbone\n'
              'SMILES: '
              '[C@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCOC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)CO, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n'
              'SMILES: CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC([*])=O, Result: '
              'True, Reason: Molecule contains the sphingosine backbone\n'
              'SMILES: '
              'CCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n'
              'SMILES: '
              'CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO[C@@H]1O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]1O)NC(=O)CCCCCCCCCCCCC\\C=C/CCCCCCCC, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n'
              'SMILES: '
              'CCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC, '
              'Result: True, Reason: Molecule contains the sphingosine '
              'backbone\n',
    'num_true_positives': 11,
    'num_false_positives': 11,
    'num_true_negatives': 0,
    'num_false_negatives': 0,
    'precision': 0.5,
    'recall': 1.0,
    'f1': 0.6666666666666666,
    'accuracy': None}