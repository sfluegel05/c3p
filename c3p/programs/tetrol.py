"""
Classifies: CHEBI:33573 tetrol
"""
from rdkit import Chem

def is_tetrol(smiles: str):
    """
    Determines if a molecule is a tetrol (a polyol that contains 4 hydroxy groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count the number of hydroxyl groups (OH)
    hydroxyl_count = 0
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'O':
            for bond in atom.GetBonds():
                neighbor = bond.GetOtherAtom(atom)
                if neighbor.GetSymbol() == 'C' and bond.GetBondType().name == 'SINGLE':
                    hydroxyl_count += 1
                    break
    
    if hydroxyl_count == 4:
        return True, "Molecule contains 4 hydroxy groups"
    else:
        return False, f"Molecule contains {hydroxyl_count} hydroxy groups"

# Example usage
smiles_list = [
    "C(\\[C@H]1[C@@H](CC([C@@H]1C/C=C\\CCCC(=O)OCC(CO)O)=O)O)=C/[C@H](CCCCC)O",
    "C([C@H]([C@@H](O)[C@H](CO)O)O)C[C@H]([C@]1(CC[C@@]2([C@]3(CC[C@@]4([C@]5(CCCC([C@@]5(CC[C@]4([C@@]3(CC[C@@]12[H])C)C)[H])(C)C)C)[H])[H])C)[H])C",
    "COC(=O)[C@]12OC[C@@]34[C@H]1[C@@H](O)C(=O)O[C@@H]3C[C@H]1C(C)=C(O)C(=O)C[C@]1(C)[C@H]4[C@@H](O)[C@@H]2O",
    "C1[C@]2([C@](/C(=C/C=C/3\\C([C@H](C[C@@H](C3)O)O)=C)/CC1)(CC[C@@]2([C@H](CC#CC(C(F)(F)F)(C(F)(F)F)O)CCCC(C([2H])([2H])[2H])(C([2H])([2H])[2H])O)[H])[H])C",
    "OCCN(CCO)c1nc(N2CCCCC2)c2nc(nc(N3CCCCC3)c2n1)N(CCO)CCO",
    "OC(C(O)C(O)C=N)C(O)C",
    "O[C@H]1[C@H](O)[C@@H](O)[C@H]2O[C@H]2[C@@H]1O",
    "COc1cc(ccc1O)[C@H](O)[C@@H](O)CO",
    "C[C@H](C[C@@H](O)[C@@H](O)C(C)(C)O)[C@H]1CC[C@@]2(C)[C@@H]3[C@@H](O)C[C@@H]4[C@]5(C[C@@]35CC[C@]12C)CCC(=O)C4(C)C",
    "C(CCCCCCCCCCCC(C(C(C(N)CO)O)O)O)CCC"
]

for smiles in smiles_list:
    result, reason = is_tetrol(smiles)
    print(f"SMILES: {smiles}\nIs Tetrol: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33573',
                          'name': 'tetrol',
                          'definition': 'A polyol that contains 4 hydroxy '
                                        'groups.',
                          'parents': ['CHEBI:26191']},
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
              'C(\\[C@H]1[C@@H](CC([C@@H]1C/C=C\\CCCC(=O)OCC(CO)O)=O)O)=C/[C@H](CCCCC)O\n'
              'Is Tetrol: False\n'
              'Reason: Molecule contains 5 hydroxy groups\n'
              '\n'
              'SMILES: '
              'C([C@H]([C@@H](O)[C@H](CO)O)O)C[C@H]([C@]1(CC[C@@]2([C@]3(CC[C@@]4([C@]5(CCCC([C@@]5(CC[C@]4([C@@]3(CC[C@@]12[H])C)C)[H])(C)C)C)[H])[H])C)[H])C\n'
              'Is Tetrol: True\n'
              'Reason: Molecule contains 4 hydroxy groups\n'
              '\n'
              'SMILES: '
              'COC(=O)[C@]12OC[C@@]34[C@H]1[C@@H](O)C(=O)O[C@@H]3C[C@H]1C(C)=C(O)C(=O)C[C@]1(C)[C@H]4[C@@H](O)[C@@H]2O\n'
              'Is Tetrol: False\n'
              'Reason: Molecule contains 7 hydroxy groups\n'
              '\n'
              'SMILES: '
              'C1[C@]2([C@](/C(=C/C=C/3\\C([C@H](C[C@@H](C3)O)O)=C)/CC1)(CC[C@@]2([C@H](CC#CC(C(F)(F)F)(C(F)(F)F)O)CCCC(C([2H])([2H])[2H])(C([2H])([2H])[2H])O)[H])[H])C\n'
              'Is Tetrol: True\n'
              'Reason: Molecule contains 4 hydroxy groups\n'
              '\n'
              'SMILES: OCCN(CCO)c1nc(N2CCCCC2)c2nc(nc(N3CCCCC3)c2n1)N(CCO)CCO\n'
              'Is Tetrol: True\n'
              'Reason: Molecule contains 4 hydroxy groups\n'
              '\n'
              'SMILES: OC(C(O)C(O)C=N)C(O)C\n'
              'Is Tetrol: True\n'
              'Reason: Molecule contains 4 hydroxy groups\n'
              '\n'
              'SMILES: O[C@H]1[C@H](O)[C@@H](O)[C@H]2O[C@H]2[C@@H]1O\n'
              'Is Tetrol: False\n'
              'Reason: Molecule contains 5 hydroxy groups\n'
              '\n'
              'SMILES: COc1cc(ccc1O)[C@H](O)[C@@H](O)CO\n'
              'Is Tetrol: False\n'
              'Reason: Molecule contains 5 hydroxy groups\n'
              '\n'
              'SMILES: '
              'C[C@H](C[C@@H](O)[C@@H](O)C(C)(C)O)[C@H]1CC[C@@]2(C)[C@@H]3[C@@H](O)C[C@@H]4[C@]5(C[C@@]35CC[C@]12C)CCC(=O)C4(C)C\n'
              'Is Tetrol: True\n'
              'Reason: Molecule contains 4 hydroxy groups\n'
              '\n'
              'SMILES: C(CCCCCCCCCCCC(C(C(C(N)CO)O)O)O)CCC\n'
              'Is Tetrol: True\n'
              'Reason: Molecule contains 4 hydroxy groups\n'
              '\n',
    'num_true_positives': 6,
    'num_false_positives': 2,
    'num_true_negatives': 8,
    'num_false_negatives': 4,
    'precision': 0.75,
    'recall': 0.6,
    'f1': 0.6666666666666665,
    'accuracy': None}