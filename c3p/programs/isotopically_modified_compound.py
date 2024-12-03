"""
Classifies: CHEBI:139358 isotopically modified compound
"""
from rdkit import Chem

def is_isotopically_modified_compound(smiles: str):
    """
    Determines if a molecule is an isotopically modified compound.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isotopically modified compound, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for isotopes in the molecule
    isotopes = []
    for atom in mol.GetAtoms():
        isotope = atom.GetIsotope()
        if isotope != 0:
            isotopes.append((atom.GetSymbol(), isotope))

    if isotopes:
        isotope_info = ", ".join([f"{symbol}-{isotope}" for symbol, isotope in isotopes])
        return True, f"Isotopically modified compound with isotopes: {isotope_info}"
    else:
        return False, "No isotopic modifications detected"

# Example usage
smiles_examples = [
    "[N+](CCO)(C)(C)C.[Cl-]",  # (11)C-choline chloride
    "O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(C([C@H](O)C(C4)([2H])[2H])([2H])[2H])[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CCC(O)=O)C)[H])[H])C",  # deoxycholic acid-2,2,4,4-d4
    "OC(=O)[C@@H](N)CCC(NC(=O)N)([2H])[2H]",  # L-citrulline-d2
    "O(C(=O)C(N([H])[H])(C(C(C([2H])([2H])[2H])(C([2H])([2H])[2H])[2H])([2H])[2H])[2H])[H]",  # leucine-d10
    "[H][C@]1([18F])C(O)O[C@H](CO)[C@@H](O)[C@@H]1O",  # 2-deoxy-2-((18)F)fluoro-D-glucopyranose
    "[C@H](OC(CCCCCCCCC)=O)(C[N+](C([2H])([2H])[2H])(C)C)CC(=O)[O-]",  # decanoyl-L-carnitine-d3
    "[N+]=1(C=CC=C(C1)C(=O)N)C([2H])([2H])[2H]",  # 1-methyl nicotinamide-d3
    "S1[C@H]([C@]2(NC(=O)N[C@]2(C1)[H])[H])CCCC(C(O)=O)([2H])[2H]",  # d2-biotin
    "O[13C](=O)[13CH2][15NH2]",  # glycine-13C2,15N
    "C[Si](C(C(C(S(O)(=O)=O)([2H])[2H])([2H])[2H])([2H])[2H])(C)C",  # 3-(trimethylsilyl)-1-propanesulfonic acid-d6
    "C(C(C(C(=O)O)([2H])[2H])([2H])[2H])(O)=O",  # succinic acid-d4
    "[13CH3]N1C=NC2=C1C(=O)N([13CH3])C(=O)N2[13CH3]",  # caffeine-(trimethyl-(13)C3)
    "C([C@@H](CC([O-])=O)OC(CCCC)=O)[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]",  # valeryl-L-carnitine-d9
    "[H][18O][H]",  # ((18)O)water
    "C(C(C([2H])([2H])[2H])([2H])N([H])[H])(=O)O[H]"  # alanine-2,3,3,3-d4
]

for smiles in smiles_examples:
    result, reason = is_isotopically_modified_compound(smiles)
    print(f"SMILES: {smiles}\nResult: {result}\nReason: {reason}\n")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:139358',
                          'name': 'isotopically modified compound',
                          'definition': 'Any molecular entity in which the '
                                        'isotopic ratio of nuclides for at '
                                        'least one element deviates measurably '
                                        'from that occurring in nature. The '
                                        'term includes both isotopically '
                                        'substituted compounds (in which '
                                        'essentially all the molecules of the '
                                        'compound have only the indicated '
                                        'nuclide(s) at each designated '
                                        'position) and isotopically labeled '
                                        'compounds (a formal mixture of an '
                                        'isotopically unmodified compound with '
                                        'one or more analogous isotopically '
                                        'substituted compound(s).',
                          'parents': ['CHEBI:23367']},
    'config': {   'llm_model_name': 'lbl/gpt-4o',
                  'accuracy_threshold': 0.95,
                  'max_attempts': 5,
                  'max_negative': 20,
                  'test_proportion': 0.1},
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '[19:22:34] SMILES Parse Error: syntax error while parsing: '
             'O[C@@H]1C[C@@H](C([O-])=O)\\[N+](C1)=C/C=C1C[C@H](NC(=C\x01)C(O)=O)C(O)=O\n'
             '[19:22:34] SMILES Parse Error: Failed parsing SMILES '
             "'O[C@@H]1C[C@@H](C([O-])=O)\\[N+](C1)=C/C=C1C[C@H](NC(=C\x01)C(O)=O)C(O)=O' "
             'for input: '
             "'O[C@@H]1C[C@@H](C([O-])=O)\\[N+](C1)=C/C=C1C[C@H](NC(=C\x01)C(O)=O)C(O)=O'\n",
    'stdout': 'SMILES: [N+](CCO)(C)(C)C.[Cl-]\n'
              'Result: False\n'
              'Reason: No isotopic modifications detected\n'
              '\n'
              'SMILES: '
              'O[C@@H]1[C@]2([C@]([C@]3([C@@]([C@@]4([C@](CC3)(C([C@H](O)C(C4)([2H])[2H])([2H])[2H])[H])C)(C1)[H])[H])(CC[C@@]2([C@@H](CCC(O)=O)C)[H])[H])C\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2, '
              'H-2, H-2\n'
              '\n'
              'SMILES: OC(=O)[C@@H](N)CCC(NC(=O)N)([2H])[2H]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2\n'
              '\n'
              'SMILES: '
              'O(C(=O)C(N([H])[H])(C(C(C([2H])([2H])[2H])(C([2H])([2H])[2H])[2H])([2H])[2H])[2H])[H]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2, '
              'H-2, H-2, H-2, H-2, H-2, H-2, H-2, H-2\n'
              '\n'
              'SMILES: [H][C@]1([18F])C(O)O[C@H](CO)[C@@H](O)[C@@H]1O\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: F-18\n'
              '\n'
              'SMILES: '
              '[C@H](OC(CCCCCCCCC)=O)(C[N+](C([2H])([2H])[2H])(C)C)CC(=O)[O-]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2, '
              'H-2\n'
              '\n'
              'SMILES: [N+]=1(C=CC=C(C1)C(=O)N)C([2H])([2H])[2H]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2, '
              'H-2\n'
              '\n'
              'SMILES: '
              'S1[C@H]([C@]2(NC(=O)N[C@]2(C1)[H])[H])CCCC(C(O)=O)([2H])[2H]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2\n'
              '\n'
              'SMILES: O[13C](=O)[13CH2][15NH2]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: C-13, '
              'C-13, N-15\n'
              '\n'
              'SMILES: '
              'C[Si](C(C(C(S(O)(=O)=O)([2H])[2H])([2H])[2H])([2H])[2H])(C)C\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2, '
              'H-2, H-2, H-2, H-2\n'
              '\n'
              'SMILES: C(C(C(C(=O)O)([2H])[2H])([2H])[2H])(O)=O\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2, '
              'H-2, H-2\n'
              '\n'
              'SMILES: [13CH3]N1C=NC2=C1C(=O)N([13CH3])C(=O)N2[13CH3]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: C-13, '
              'C-13, C-13\n'
              '\n'
              'SMILES: '
              'C([C@@H](CC([O-])=O)OC(CCCC)=O)[N+](C([2H])([2H])[2H])(C([2H])([2H])[2H])C([2H])([2H])[2H]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2, '
              'H-2, H-2, H-2, H-2, H-2, H-2, H-2\n'
              '\n'
              'SMILES: [H][18O][H]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: O-18\n'
              '\n'
              'SMILES: C(C(C([2H])([2H])[2H])([2H])N([H])[H])(=O)O[H]\n'
              'Result: True\n'
              'Reason: Isotopically modified compound with isotopes: H-2, H-2, '
              'H-2, H-2\n'
              '\n',
    'num_true_positives': 14,
    'num_false_positives': 0,
    'num_true_negatives': 15,
    'num_false_negatives': 1,
    'precision': 1.0,
    'recall': 0.9333333333333333,
    'f1': 0.9655172413793104,
    'accuracy': None}