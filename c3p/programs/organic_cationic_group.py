"""
Classifies: CHEBI:64769 organic cationic group
"""
from rdkit import Chem

def is_organic_cationic_group(smiles: str):
    """
    Determines if a molecule is an organic cationic group (a cationic group that contains carbon).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organic cationic group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for cationic charge
    if not any(atom.GetFormalCharge() > 0 for atom in mol.GetAtoms()):
        return False, "No cationic charge found"

    # Check for the presence of carbon atoms
    if not any(atom.GetSymbol() == 'C' for atom in mol.GetAtoms()):
        return False, "No carbon atoms found"

    return True, "Organic cationic group found"

# Example usage
smiles_list = [
    "CC(C)C[C@H]([NH3+])C(=O)N[C@@H]([*])C(-*)=O",
    "C(*)(=O)[C@@H](N*)CC=1[NH+]=CNC1",
    "C(=O)([C@@H]([NH3+])CCSC)N[C@H](C(=O)*)C(C)C",
    "[C@@H](C(*)=O)(CC(CC[NH3+])O)N*",
    "O=C(*)[C@@H]([NH3+])C(C)C",
    "[NH+](=C/C=C(/C=C/C=C(/C=C/C1=C(CCCC1(C)C)C)\C)\C)/CCCC[C@@H](C(*)=O)N*",
    "C(*)(=O)[C@@H]([NH3+])CC1=CC=CC=C1",
    "C([C@H](CO)[N+](C)(C)C)(=O)N1[C@H](C(=O)N[C@H](C(=O)*)CCCC[NH3+])CCC1",
    "[C@@H](C(*)=O)(CC(C(C[NH3+])O)O)N*",
    "C([C@H](CO)[NH3+])(=O)*",
    "[NH3+]CC(-*)=O",
    "O=C(*)[C@@H](N*)CCCNC(=[NH2+])NC(C(C)=O)O",
    "O=C(*)[C@@H]([NH3+])[C@H](O)C",
    "C(*)(=O)[C@@H]([NH2+]*)CC[S@](=O)C"
]

for smiles in smiles_list:
    is_cationic, reason = is_organic_cationic_group(smiles)
    print(f"SMILES: {smiles} -> {is_cationic}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:64769',
                          'name': 'organic cationic group',
                          'definition': 'A cationic group that contains '
                                        'carbon.',
                          'parents': ['CHEBI:64766']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: CC(C)C[C@H]([NH3+])C(=O)N[C@@H]([*])C(-*)=O -> True, '
              'Organic cationic group found\n'
              'SMILES: C(*)(=O)[C@@H](N*)CC=1[NH+]=CNC1 -> True, Organic '
              'cationic group found\n'
              'SMILES: C(=O)([C@@H]([NH3+])CCSC)N[C@H](C(=O)*)C(C)C -> True, '
              'Organic cationic group found\n'
              'SMILES: [C@@H](C(*)=O)(CC(CC[NH3+])O)N* -> True, Organic '
              'cationic group found\n'
              'SMILES: O=C(*)[C@@H]([NH3+])C(C)C -> True, Organic cationic '
              'group found\n'
              'SMILES: '
              '[NH+](=C/C=C(/C=C/C=C(/C=C/C1=C(CCCC1(C)C)C)\\C)\\C)/CCCC[C@@H](C(*)=O)N* '
              '-> True, Organic cationic group found\n'
              'SMILES: C(*)(=O)[C@@H]([NH3+])CC1=CC=CC=C1 -> True, Organic '
              'cationic group found\n'
              'SMILES: '
              'C([C@H](CO)[N+](C)(C)C)(=O)N1[C@H](C(=O)N[C@H](C(=O)*)CCCC[NH3+])CCC1 '
              '-> True, Organic cationic group found\n'
              'SMILES: [C@@H](C(*)=O)(CC(C(C[NH3+])O)O)N* -> True, Organic '
              'cationic group found\n'
              'SMILES: C([C@H](CO)[NH3+])(=O)* -> True, Organic cationic group '
              'found\n'
              'SMILES: [NH3+]CC(-*)=O -> True, Organic cationic group found\n'
              'SMILES: O=C(*)[C@@H](N*)CCCNC(=[NH2+])NC(C(C)=O)O -> True, '
              'Organic cationic group found\n'
              'SMILES: O=C(*)[C@@H]([NH3+])[C@H](O)C -> True, Organic cationic '
              'group found\n'
              'SMILES: C(*)(=O)[C@@H]([NH2+]*)CC[S@](=O)C -> True, Organic '
              'cationic group found\n',
    'num_true_positives': 14,
    'num_false_positives': 1,
    'num_true_negatives': 13,
    'num_false_negatives': 0,
    'precision': 0.9333333333333333,
    'recall': 1.0,
    'f1': 0.9655172413793104,
    'accuracy': None}