"""
Classifies: CHEBI:27207 univalent carboacyl group
"""
from rdkit import Chem

def is_univalent_carboacyl_group(smiles: str):
    """
    Determines if a molecule is a univalent carboacyl group (formed by loss of OH from the carboxy group of a carboxylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a univalent carboacyl group, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule contains a carbonyl group (C=O)
    carbonyl = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'C':
            for neighbor in atom.GetNeighbors():
                if neighbor.GetSymbol() == 'O' and mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                    carbonyl = True
                    break
            if carbonyl:
                break

    if not carbonyl:
        return False, "No carbonyl group (C=O) found"

    # Check if the molecule has a valence of 1 (univalent) and contains a '*' atom
    univalent = False
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == '*':
            univalent = True
            break

    if not univalent:
        return False, "Molecule is not univalent (missing * atom)"

    return True, "Molecule is a univalent carboacyl group"

# Test examples
smiles_list = [
    "CCCCC\\C=C/C\\C=C/CC#CCCCCC(-*)=O",
    "C1=C(C(=C(C=C1)Cl)C(=O)*)Cl",
    "[*]C(=O)CC(-*)=O",
    "C=1(/C=C/C(=O)*)C=C(OC)C(O)=C(OC)C1",
    "C1[C@@H](SSC1)CCCCC(*)=O",
    "C(C(*)=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC",
    "C(=O)(*)CC1=CC(=C(C(=C1)Br)O)[N+]([O-])=O",
    "C(*)(=O)[C@H](N)C",
    "CCC\\C=C\\C\\C=C/CCCCCCCC(-*)=O",
    "C(CCCCCCC#CCCCCCCCC(=O)*)C",
    "C1(=CNC2=C1C=CC=C2)C[C@@H](C(=O)*)N",
    "C(C1=CC(=C(C(=C1)I)O)[N+]([O-])=O)C(=O)*",
    "[NH3+][C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)*)CC(=O)[O-])CC(=O)[O-])CC(=O)[O-]",
    "[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O[C@H](C)C(=O)*)N)O",
    "CC\\C=C/C\\C=C/C\\C=C/CC#CCCCCC(-*)=O",
    "C([C@H](C(*)=O)N)C(N)=O",
    "C([C@]([C@]([C@]([C@](CO)([H])O)([H])O)(O[C@H](C)C(=O)*)[H])([H])N)=O",
    "OC(=O)[C@H](N)CCC(=O)*",
    "CCCCCCC(=C(CCCCCCCCCC(=O)*)[H])[H]",
    "C([C@H](C(*)=O)N)[Se]C",
    "[C@H](NC(=O)/C(/C1=CSC(=N1)N)=N\\OC(C(=O)O)(C)C)([C@H](C)NS(O)(=O)=O)C(=O)*",
    "C(C(*)=O)CCCCC",
    "C(CCCCCBr)(*)=O",
    "C(C=CCCCCC(=O)*)CCCCCCCCCC",
    "[C@@H]1(O[C@H]([C@@H]([C@@H]1O)O)N2C=CC(NC2=O)=O)COP(*)(O)=O",
    "C(*)(C(=O)O)=O",
    "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCC(-*)=O"
]

for smiles in smiles_list:
    result, reason = is_univalent_carboacyl_group(smiles)
    print(f"SMILES: {smiles} -> {result}, {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27207',
                          'name': 'univalent carboacyl group',
                          'definition': 'A univalent carboacyl group is a '
                                        'group formed by loss of OH from the '
                                        'carboxy group of a carboxylic acid.',
                          'parents': ['CHEBI:37838']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 2,
    'success': True,
    'best': True,
    'error': '[21:15:09] WARNING: not removing hydrogen atom without '
             'neighbors\n'
             '[21:15:09] WARNING: not removing hydrogen atom without '
             'neighbors\n',
    'stdout': 'SMILES: CCCCC\\C=C/C\\C=C/CC#CCCCCC(-*)=O -> True, Molecule is '
              'a univalent carboacyl group\n'
              'SMILES: C1=C(C(=C(C=C1)Cl)C(=O)*)Cl -> True, Molecule is a '
              'univalent carboacyl group\n'
              'SMILES: [*]C(=O)CC(-*)=O -> True, Molecule is a univalent '
              'carboacyl group\n'
              'SMILES: C=1(/C=C/C(=O)*)C=C(OC)C(O)=C(OC)C1 -> True, Molecule '
              'is a univalent carboacyl group\n'
              'SMILES: C1[C@@H](SSC1)CCCCC(*)=O -> True, Molecule is a '
              'univalent carboacyl group\n'
              'SMILES: C(C(*)=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC -> True, '
              'Molecule is a univalent carboacyl group\n'
              'SMILES: C(=O)(*)CC1=CC(=C(C(=C1)Br)O)[N+]([O-])=O -> True, '
              'Molecule is a univalent carboacyl group\n'
              'SMILES: C(*)(=O)[C@H](N)C -> True, Molecule is a univalent '
              'carboacyl group\n'
              'SMILES: CCC\\C=C\\C\\C=C/CCCCCCCC(-*)=O -> True, Molecule is a '
              'univalent carboacyl group\n'
              'SMILES: C(CCCCCCC#CCCCCCCCC(=O)*)C -> True, Molecule is a '
              'univalent carboacyl group\n'
              'SMILES: C1(=CNC2=C1C=CC=C2)C[C@@H](C(=O)*)N -> True, Molecule '
              'is a univalent carboacyl group\n'
              'SMILES: C(C1=CC(=C(C(=C1)I)O)[N+]([O-])=O)C(=O)* -> True, '
              'Molecule is a univalent carboacyl group\n'
              'SMILES: '
              '[NH3+][C@H](C(=O)N[C@H](C(=O)N[C@H](C(=O)*)CC(=O)[O-])CC(=O)[O-])CC(=O)[O-] '
              '-> True, Molecule is a univalent carboacyl group\n'
              'SMILES: '
              '[C@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O[C@H](C)C(=O)*)N)O -> '
              'True, Molecule is a univalent carboacyl group\n'
              'SMILES: CC\\C=C/C\\C=C/C\\C=C/CC#CCCCCC(-*)=O -> True, Molecule '
              'is a univalent carboacyl group\n'
              'SMILES: C([C@H](C(*)=O)N)C(N)=O -> True, Molecule is a '
              'univalent carboacyl group\n'
              'SMILES: '
              'C([C@]([C@]([C@]([C@](CO)([H])O)([H])O)(O[C@H](C)C(=O)*)[H])([H])N)=O '
              '-> True, Molecule is a univalent carboacyl group\n'
              'SMILES: OC(=O)[C@H](N)CCC(=O)* -> True, Molecule is a univalent '
              'carboacyl group\n'
              'SMILES: CCCCCCC(=C(CCCCCCCCCC(=O)*)[H])[H] -> True, Molecule is '
              'a univalent carboacyl group\n'
              'SMILES: C([C@H](C(*)=O)N)[Se]C -> True, Molecule is a univalent '
              'carboacyl group\n'
              'SMILES: '
              '[C@H](NC(=O)/C(/C1=CSC(=N1)N)=N\\OC(C(=O)O)(C)C)([C@H](C)NS(O)(=O)=O)C(=O)* '
              '-> True, Molecule is a univalent carboacyl group\n'
              'SMILES: C(C(*)=O)CCCCC -> True, Molecule is a univalent '
              'carboacyl group\n'
              'SMILES: C(CCCCCBr)(*)=O -> True, Molecule is a univalent '
              'carboacyl group\n'
              'SMILES: C(C=CCCCCC(=O)*)CCCCCCCCCC -> True, Molecule is a '
              'univalent carboacyl group\n'
              'SMILES: '
              '[C@@H]1(O[C@H]([C@@H]([C@@H]1O)O)N2C=CC(NC2=O)=O)COP(*)(O)=O -> '
              'True, Molecule is a univalent carboacyl group\n'
              'SMILES: C(*)(C(=O)O)=O -> True, Molecule is a univalent '
              'carboacyl group\n'
              'SMILES: CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCC(-*)=O -> True, '
              'Molecule is a univalent carboacyl group\n',
    'num_true_positives': 27,
    'num_false_positives': 1,
    'num_true_negatives': 19,
    'num_false_negatives': 0,
    'precision': 0.9642857142857143,
    'recall': 1.0,
    'f1': 0.9818181818181818,
    'accuracy': None}