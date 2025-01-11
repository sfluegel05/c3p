"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
from rdkit import Chem

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule is a 3beta-hydroxy-Delta(5)-steroid based on its SMILES string.
    A 3beta-hydroxy-Delta(5)-steroid contains a 3beta-hydroxy group and a double bond between positions 5 and 6.

    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3beta-hydroxy-Delta(5)-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # General SMARTS pattern for 3beta-hydroxy group
    hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)[C;R1]")  # A 3beta-hydroxy group on a cyclic structure
    
    # General SMARTS pattern for Delta(5) double bond in steroid
    delta5_pattern = Chem.MolFromSmarts("C=C[C;R]")  # A double bond in a cyclic structure

    # Check for 3beta-hydroxy group
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3beta-hydroxy group found"
    
    # Check for Delta(5) double bond
    # We manually check by identifying positions between rings where double bonds occur
    carbon_bond_info = [(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
                        for bond in mol.GetBonds() if bond.GetBondTypeAsDouble() == 2.0]
    delta_5_present = any((5 in pair and 6 in pair) or (pair == sorted(pair)[::-1] and all(5 <= idx < 7 for idx in pair)) for pair in carbon_bond_info)

    if not delta_5_present:
        return False, "No Delta(5) double bond found (C5-C6 double bond)"

    return True, "Contains 3beta-hydroxy group with Delta(5) double bond (C5-C6 double bond)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:1722',
                          'name': '3beta-hydroxy-Delta(5)-steroid',
                          'definition': 'Any 3beta-hydroxy-steroid that '
                                        'contains a double bond between '
                                        'positions 5 and 6.',
                          'parents': ['CHEBI:36836'],
                          'xrefs': [   'KEGG:C03836',
                                       'MetaCyc:3b-hydroxy-D5-steroids'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@H](O)C(C)C '
               'NAME: (24S)-7alpha,24-dihydroxycholesterol REASON: MISSED No '
               '3beta-hydroxy group found\n'
               ' * SMILES: '
               '[H][C@]12C[C@@]3([H])[C@]4([H])CC=C5C[C@@H](O)CC[C@]5(C)[C@@]4([H])CC[C@]3(C)[C@@]1(O)[C@H](C)[C@@]1(CC[C@@H](C)CO1)O2 '
               'NAME: pennogenin REASON: MISSED No 3beta-hydroxy group found\n'
               ' * SMILES: '
               'C[C@H](CCC[C@@H](C)C(O)=O)[C@H]1CC[C@H]2[C@@H]3CC=C4[C@@H](O)[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C '
               'NAME: (25R)-3beta,4beta-dihydroxycholest-5-en-26-oic acid '
               'REASON: MISSED No 3beta-hydroxy group found\n'
               ' * SMILES: '
               '[H][C@@]12CC[C@H](C(C)=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])[C@@H](O)C=C2C[C@@H](O)CC[C@]12C '
               'NAME: 7beta-hydroxypregnenolone REASON: MISSED No '
               '3beta-hydroxy group found\n'
               ' * SMILES: '
               'O([C@@]1(C(C=2[C@@]([C@@]3(C([C@]4([C@@]([C@](C(C4([H])[H])([H])[H])([C@](C([H])([H])[H])(/C(=C(/[C@](C(C([H])([H])[H])(C([H])([H])[H])[H])(C(C([H])([H])[H])([H])[H])[H])\\[H])/[H])[H])[H])(C(C3([H])[H])([H])[H])C([H])([H])[H])[H])=C(C2[H])[H])[H])(C(C1([H])[H])([H])[H])C([H])([H])[H])([H])[H])[H])[H] '
               'NAME: delta7-stigmasterol REASON: MISSED No 3beta-hydroxy '
               'group found\n'
               ' * SMILES: '
               'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C(C=C2C[C@H](C1)O)=O)[H])(CC[C@@]4([C@H](C)CCC[C@@H](C)CO)[H])[H])C)[H])C '
               'NAME: (25R)-3beta,26-dihydroxycholest-5-en-7-one REASON: '
               'MISSED No 3beta-hydroxy group found\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCC(O)C(C)C '
               'NAME: 7alpha,24-dihydroxycholesterol REASON: MISSED No '
               '3beta-hydroxy group found\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCCC(C)CO '
               'NAME: 7alpha,26-dihydroxycholesterol REASON: MISSED No '
               '3beta-hydroxy group found\n'
               ' * SMILES: '
               '[H][C@@]12CC=C3C[C@@H](O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)C(=O)CC[C@@]21[H] '
               'NAME: dehydroepiandrosterone REASON: MISSED No 3beta-hydroxy '
               'group found\n'
               ' * SMILES: '
               '[C@@]12([C@]3([C@](CC[C@@]1([C@@]4(C(C[C@H](CC4)O)=CC2)C)[H])([C@](CC3)([C@H](C)CC[C@@H](C(C)C)OS(O)(=O)=O)[H])C)[H])[H] '
               'NAME: (24S)-hydroxycholesterol 24-sulfate REASON: MISSED No '
               '3beta-hydroxy group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@@H](C(=O)[C@@]2([C@]1(C(=C)C[C@H]3[C@]4([C@@H](C(C(=O)CC4)(C)C)[C@H](C[C@]23C)O)C)C)C(=O)OC)C',
                                     'name': 'Terreustoxin L',
                                     'reason': 'No Delta(5) double bond found '
                                               '(C5-C6 double bond)'},
                                 {   'smiles': '[H][C@]12CCC[C@H](C)[C@@]1(C)C[C@@H](CC2)C(C)C',
                                     'name': 'eremophilane',
                                     'reason': 'No 3beta-hydroxy group found'},
                                 {   'smiles': 'N.C[C@]12CC[C@@](C)(C[C@H]1C1=CC(=O)[C@@H]3[C@@]4(C)CC[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)C(O)=O)C(O)=O)C(C)(C)[C@@H]4CC[C@@]3(C)[C@]1(C)CC2)C(O)=O',
                                     'name': 'Monoammonium glycyrrhizinate',
                                     'reason': 'No Delta(5) double bond found '
                                               '(C5-C6 double bond)'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@H]2[C@H](O[C@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O)O[C@@H]1CO)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)[C@@H]4O[C@@H]6OC[C@@H](O)[C@H](O)[C@H]6O)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[(2R,3S,4R,5R,6R)-5-Acetamido-6-hydroxy-2-(hydroxymethyl)-4-[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-4-hydroxy-5-[(2S,3S,4S,5R,6R)-5-hydroxy-6-(hydroxymethyl)-4-[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3-[(2S,3R,4S,5R)-3,4,5-trihydroxyoxan-2-yl]oxyoxan-2-yl]oxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No Delta(5) double bond found '
                                               '(C5-C6 double bond)'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1NC=NC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC2=CC=CC=C2)CCC(=O)N',
                                     'name': 'Phe-Gln-His',
                                     'reason': 'No 3beta-hydroxy group found'},
                                 {   'smiles': 'CN(C)CC(=O)N[C@H]1CC[C@@H](O[C@H]1CO)CCN2C=C(N=N2)C3=CC(=CC=C3)OC',
                                     'name': '2-(dimethylamino)-N-[(2R,3S,6R)-2-(hydroxymethyl)-6-[2-[4-(3-methoxyphenyl)-1-triazolyl]ethyl]-3-oxanyl]acetamide',
                                     'reason': 'No Delta(5) double bond found '
                                               '(C5-C6 double bond)'},
                                 {   'smiles': 'CC1=CC=C(C=C1)CS(=O)(=O)C2=NC=C(C(=N2)C(=O)NC3=NN=C(S3)C(C)C)Cl',
                                     'name': '5-chloro-2-[(4-methylphenyl)methylsulfonyl]-N-(5-propan-2-yl-1,3,4-thiadiazol-2-yl)-4-pyrimidinecarboxamide',
                                     'reason': 'No 3beta-hydroxy group found'},
                                 {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\CCCCC)([O-])=O',
                                     'name': 'PC(15:1(9Z)/15:1(9Z))',
                                     'reason': 'No 3beta-hydroxy group found'},
                                 {   'smiles': 'OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)C)[C@H](CC)C',
                                     'name': 'Asp-Ala-Ile',
                                     'reason': 'No 3beta-hydroxy group found'},
                                 {   'smiles': 'O=C(OC[C@]1(O)[C@@H]2C=C3C=C[C@@H]4[C@@]([C@]3(C2)CC1)(CC[C@H]([C@]4(CO)C)O)C)C',
                                     'name': 'Aphidicolin A42',
                                     'reason': 'No Delta(5) double bond found '
                                               '(C5-C6 double bond)'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@H](O)C(C)C',
                                      'name': '(24S)-7alpha,24-dihydroxycholesterol',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': 'C[C@H](CCC[C@@H](C)C(O)=O)[C@H]1CC[C@H]2[C@@H]3CC=C4[C@@H](O)[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C',
                                      'name': '(25R)-3beta,4beta-dihydroxycholest-5-en-26-oic '
                                              'acid',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': '[H][C@@]12CC[C@H](C(C)=O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])[C@@H](O)C=C2C[C@@H](O)CC[C@]12C',
                                      'name': '7beta-hydroxypregnenolone',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': 'O([C@@]1(C(C=2[C@@]([C@@]3(C([C@]4([C@@]([C@](C(C4([H])[H])([H])[H])([C@](C([H])([H])[H])(/C(=C(/[C@](C(C([H])([H])[H])(C([H])([H])[H])[H])(C(C([H])([H])[H])([H])[H])[H])\\[H])/[H])[H])[H])(C(C3([H])[H])([H])[H])C([H])([H])[H])[H])=C(C2[H])[H])[H])(C(C1([H])[H])([H])[H])C([H])([H])[H])([H])[H])[H])[H]',
                                      'name': 'delta7-stigmasterol',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C(C=C2C[C@H](C1)O)=O)[H])(CC[C@@]4([C@H](C)CCC[C@@H](C)CO)[H])[H])C)[H])C',
                                      'name': '(25R)-3beta,26-dihydroxycholest-5-en-7-one',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCC(O)C(C)C',
                                      'name': '7alpha,24-dihydroxycholesterol',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCCC(C)CO',
                                      'name': '7alpha,26-dihydroxycholesterol',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': '[H][C@@]12CC=C3C[C@@H](O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)C(=O)CC[C@@]21[H]',
                                      'name': 'dehydroepiandrosterone',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': '[C@@]12([C@]3([C@](CC[C@@]1([C@@]4(C(C[C@H](CC4)O)=CC2)C)[H])([C@](CC3)([C@H](C)CC[C@@H](C(C)C)OS(O)(=O)=O)[H])C)[H])[H]',
                                      'name': '(24S)-hydroxycholesterol '
                                              '24-sulfate',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'},
                                  {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C(C=C2C[C@H](C1)O)=O)[H])(CC[C@@]4([C@H](C)CCC[C@@H](C)C([H])=O)[H])[H])C)[H])C',
                                      'name': '(25R)-3beta-hydroxycholest-5-en-7-one-26-al',
                                      'reason': 'No Delta(5) double bond found '
                                                '(C5-C6 double bond)'}],
    'attempt': 2,
    'reasoning': 'From the outcomes, it seems the previous code did not '
                 'correctly identify the 3beta-hydroxy group or the Delta(5) '
                 'double bond. This could be due to overly specific SMARTS '
                 'patterns that fail to capture the necessary chemical '
                 'features in the provided steroid examples.\n'
                 '\n'
                 'To improve the classification, the SMARTS patterns must be '
                 'designed to accurately identify:\n'
                 '\n'
                 '1. **3beta-hydroxy group**: The program must ensure the '
                 '3beta-hydroxy group, which is a hydroxyl group attached at '
                 'the beta orientation to the third carbon of the steroid '
                 'structure.\n'
                 '   \n'
                 '2. **Delta(5) double bond**: This is a double bond between '
                 'the 5th and 6th carbon atoms in the steroid, a common '
                 'feature in Delta(5)-steroids.\n'
                 '\n'
                 'To improve the accuracy, a more generalized approach using '
                 'broader SMARTS patterns that effectively capture these '
                 "structures should be employed. Here's the refined code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 5,
    'num_false_positives': 100,
    'num_true_negatives': 3842,
    'num_false_negatives': 87,
    'num_negatives': None,
    'precision': 0.047619047619047616,
    'recall': 0.05434782608695652,
    'f1': 0.05076142131979696,
    'accuracy': 0.9536440257808627,
    'negative_predictive_value': 0.9778569610587936}