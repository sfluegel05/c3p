"""
Classifies: CHEBI:51702 enoate ester
"""
from rdkit import Chem

def is_enoate_ester(smiles: str):
    """
    Determines if a molecule is an enoate ester (alpha,beta-unsaturated carboxylic ester).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enoate ester, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ester functional group
    ester_found = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C' and atom.GetDegree() == 3:
            neighbors = atom.GetNeighbors()
            oxygen_count = sum(1 for n in neighbors if n.GetSymbol() == 'O')
            carbon_count = sum(1 for n in neighbors if n.GetSymbol() == 'C')
            if oxygen_count == 2 and carbon_count == 1:
                ester_found = True
                ester_carbon = atom
                break

    if not ester_found:
        return False, "No ester functional group found"

    # Check for alpha,beta-unsaturated carbonyl system
    alpha_beta_unsaturated_found = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if (begin_atom.GetSymbol() == 'C' and end_atom.GetSymbol() == 'C'):
                begin_neighbors = [n for n in begin_atom.GetNeighbors() if n != end_atom]
                end_neighbors = [n for n in end_atom.GetNeighbors() if n != begin_atom]
                if any(n.GetSymbol() == 'C' and n.GetIdx() == ester_carbon.GetIdx() for n in begin_neighbors + end_neighbors):
                    if any(n.GetSymbol() == 'C' and n.GetIdx() != ester_carbon.GetIdx() for n in begin_neighbors + end_neighbors):
                        alpha_beta_unsaturated_found = True
                        break

    if not alpha_beta_unsaturated_found:
        return False, "No alpha,beta-unsaturated carbonyl system found"

    return True, "Molecule is an enoate ester"

# Examples
smiles_list = [
    "O(C(=O)C=1C(C(=C(NC1C2=CC=CC=C2)C)C(OCC)=O)C#CC3=CC=CC=C3)CC4=CC=CC=C4",
    "C[C@H]1C[C@@]2(O)[C@H]([C@H]1OC(=O)\\C=C\\c1ccccc1)[C@@H](O)\\C(C)=C/CC[C@H]1[C@@H](\\C=C(C)\\C2=O)C1(C)C",
    "O=C1O[C@@H]([C@@H]2O[C@@H]2C=CC(C=C1)=O)C",
    "CCCC(CC)COC(\\C(=C(\\C1=CC=C(C=C1)OC)/[H])\\[H])=O",
    "O(C(CCC=C(C)C)(C)C=C)C(=O)/C=C/C1=CC=CC=C1"
]

for smiles in smiles_list:
    result, reason = is_enoate_ester(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51702',
                          'name': 'enoate ester',
                          'definition': 'An alpha,beta-unsaturated carboxylic '
                                        'ester of general formula '
                                        'R(1)R(2)C=CR(3)-C(=O)OR(4) (R(4) =/= '
                                        'H) in which the ester C=O function is '
                                        'conjugated to a C=C double bond at '
                                        'the alpha,beta position.',
                          'parents': ['CHEBI:51737', 'CHEBI:78840']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': 'SMILES: '
              'O(C(=O)C=1C(C(=C(NC1C2=CC=CC=C2)C)C(OCC)=O)C#CC3=CC=CC=C3)CC4=CC=CC=C4 '
              '-> True, Molecule is an enoate ester\n'
              'SMILES: '
              'C[C@H]1C[C@@]2(O)[C@H]([C@H]1OC(=O)\\C=C\\c1ccccc1)[C@@H](O)\\C(C)=C/CC[C@H]1[C@@H](\\C=C(C)\\C2=O)C1(C)C '
              '-> True, Molecule is an enoate ester\n'
              'SMILES: O=C1O[C@@H]([C@@H]2O[C@@H]2C=CC(C=C1)=O)C -> True, '
              'Molecule is an enoate ester\n'
              'SMILES: CCCC(CC)COC(\\C(=C(\\C1=CC=C(C=C1)OC)/[H])\\[H])=O -> '
              'True, Molecule is an enoate ester\n'
              'SMILES: O(C(CCC=C(C)C)(C)C=C)C(=O)/C=C/C1=CC=CC=C1 -> True, '
              'Molecule is an enoate ester\n',
    'num_true_positives': 33,
    'num_false_positives': 0,
    'num_true_negatives': 20,
    'num_false_negatives': 11,
    'precision': 1.0,
    'recall': 0.75,
    'f1': 0.8571428571428571,
    'accuracy': None}