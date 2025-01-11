"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    Long-chain fatty alcohols contain a long hydrocarbon chain (typically 14 or more carbons)
    with an alcohol (-OH) group at the end or within the chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for alcohol group (oxygen connected to a carbon, checking for both primary, secondary)
    alcohol_pattern_primary = Chem.MolFromSmarts('[CX4][OX2H]')
    alcohol_pattern_secondary = Chem.MolFromSmarts('[CX3][OX2H]')
    alcohol_pattern_tertiary = Chem.MolFromSmarts('[CX2][OX2H]')
    if not (mol.HasSubstructMatch(alcohol_pattern_primary) or 
            mol.HasSubstructMatch(alcohol_pattern_secondary) or 
            mol.HasSubstructMatch(alcohol_pattern_tertiary)):
        return False, "No proper alcohol group found"

    # Check for total carbon count (at least 14)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 14:
        return False, "Not enough total carbons for a long hydrocarbon chain"

    # Success: Detected a long-chain fatty alcohol
    return True, "Contains long hydrocarbon chain with an alcohol group"

# Test the function with a sample SMILES string from the list:
print(is_long_chain_fatty_alcohol('CCCCCCCCCC/C=C/CCCCCCCCCCCCO'))  # Example: (13Z)-docosen-1-ol


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17135',
                          'name': 'long-chain fatty alcohol',
                          'definition': 'A fatty alcohol with a chain length '
                                        'ranging from C13 to C22.',
                          'parents': ['CHEBI:24026'],
                          'xrefs': [   'KEGG:C00339',
                                       'MetaCyc:Long-chain-alcohols',
                                       'PMID:19155068',
                                       'PMID:22436224',
                                       'PMID:22444994'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.2857142857142857 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: CCCCCCCC\\C=C/CCCCCCCCCCCCO NAME: '
               '(13Z)-docosen-1-ol REASON: CORRECT Contains long hydrocarbon '
               'chain with an alcohol group\n'
               ' * SMILES: CCCCCCCCCCCCCCO NAME: tetradecan-1-ol REASON: '
               'CORRECT Contains long hydrocarbon chain with an alcohol group\n'
               ' * SMILES: OCCCCCCCCCCCCC#CCCC NAME: 13-Heptadecyn-1-ol '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: CCCCCCCCC(O)CCCCCO NAME: 1,6-tetradecanediol '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: '
               'O=C(OC(C(O)C(O)C(=O)O)CO)C(C(=O)/C(=C/C(CC(CCCCC(O)CCCCCCCCCCCCC)C)C)/C)C '
               'NAME: Haplofungin E REASON: CORRECT Contains long hydrocarbon '
               'chain with an alcohol group\n'
               ' * SMILES: CC(C)CCCCCCCCCCCO NAME: 12-methyltridecan-1-ol '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: '
               'O=C(OC(C(O)C(O)C(=O)O)CO)C(C(=O)/C(=C/C(CC(CCCC(O)C(O)CCCCCCCCCCC)C)C)/C)C '
               'NAME: Haplofungin H REASON: CORRECT Contains long hydrocarbon '
               'chain with an alcohol group\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCCCCCCO NAME: docosan-1-ol REASON: '
               'CORRECT Contains long hydrocarbon chain with an alcohol group\n'
               ' * SMILES: CCCCCCCCCC(O)CCCCO NAME: 1,5-tetradecanediol '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: CCCCCCCCCCC(O)CCCO NAME: 1,4-tetradecanediol '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: CC(C)CCCCCCCCCCCCO NAME: 13-methyltetradecan-1-ol '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: '
               'O(C(=O)CCCCCCCCCCCCCC)C[C@@H](OC(=O)CC/C=C\\C/C=C\\C/C=C\\C=C\\C(O)C/C=C\\C/C=C\\CC)CO '
               'NAME: DG(15:0/22:6(4Z,7Z,10Z,12E,16Z,19Z)-OH(14)/0:0) REASON: '
               'CORRECT Contains long hydrocarbon chain with an alcohol group\n'
               ' * SMILES: CCC(C)CCCCCCCCCCCCCCCCCO NAME: 18-methylicosan-1-ol '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCCCCO NAME: icosan-1-ol REASON: '
               'CORRECT Contains long hydrocarbon chain with an alcohol group\n'
               ' * SMILES: '
               'O(C(=O)CCCCCCCCCCCCCCC)[C@H](COC(=O)CCC/C=C\\[C@H](O)/C=C/C=C\\C/C=C\\C=C\\[C@H](O)C/C=C\\CC)CO '
               'NAME: DG(22:6(5Z,8E,10Z,13Z,15E,19Z)-2OH(7S, 17S)/16:0/0:0) '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: CC(C)CCCCCCCCCCCCCCCCCO NAME: '
               '18-methylnonadecan-1-ol REASON: CORRECT Contains long '
               'hydrocarbon chain with an alcohol group\n'
               ' * SMILES: OCCCCCCCCCCCCCC#CCC=C NAME: 17-OCTADECEN-14-YN-1-OL '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: CCC(C)CCCCCCCCCCCCCO NAME: 14-methylhexadecan-1-ol '
               'REASON: CORRECT Contains long hydrocarbon chain with an '
               'alcohol group\n'
               ' * SMILES: '
               'O=C(OC(C(O)C(O)C(=O)O)CO)C(C(=O)/C(=C/C(CC(CCCCC(O)CCCCCCCCCCC)C)C)/C)C '
               'NAME: Haplofungin A REASON: CORRECT Contains long hydrocarbon '
               'chain with an alcohol group\n'
               ' * SMILES: '
               'O=C(OC(C(O)C(O)C(=O)O)CO)C(C(=O)/C(=C/C(CC(CCCCC(O)C(O)CCCCCCCCCC)C)C)/C)C '
               'NAME: Haplofungin F REASON: CORRECT Contains long hydrocarbon '
               'chain with an alcohol group\n'
               'False positives: NONE\n'
               'False negatives: SMILES: OCCCCCCCCCC/C=C/CC NAME: '
               '11E-Tetradecen-1-ol REASON: MISSED No long hydrocarbon chain '
               'found\n'
               ' * SMILES: OCCCCCC/C=C\\CCCCCCCC NAME: 7Z-hexadecen-1-ol '
               'REASON: MISSED No long hydrocarbon chain found\n'
               ' * SMILES: '
               'O=C1C(=C(CC(O)CC(O)/C=C/C=C/C(=C/C=C/C(C(O)C(CCC2=CC=CC=C2)C)C)/C)OC1(O)C)C '
               'NAME: Hyafurone A2 REASON: MISSED No proper alcohol group '
               'found\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\CO '
               'NAME: all-trans-undecaprenol REASON: MISSED No long '
               'hydrocarbon chain found\n'
               ' * SMILES: OC(CCCC/C=C\\CCCCCCCC1=C(C(O)=CC=C1)C(O)=O)C NAME: '
               '6-((Z)-14-hydroxypentadec-8-en-1-yl)salicylic acid REASON: '
               'MISSED No proper alcohol group found\n'
               ' * SMILES: '
               'O=C(C[C@@H](O)/C=C/C=C/C=C\\[C@@H]([C@H](O)[C@H]([C@@H](O)C[C@H](O)C)C)C)CC '
               'NAME: Chivotriene REASON: MISSED No proper alcohol group '
               'found\n'
               ' * SMILES: OC(CC(=O)CCCCCCCCCCCCCCCC)CCC NAME: '
               '4-Hydroxy-6-docosanone REASON: MISSED No proper alcohol group '
               'found\n'
               ' * SMILES: OC(/C(=C\\CCCCCCCCCCCCC)/C(OC)=O)C(=O)C NAME: '
               'Methyl (2E)-2-(1-hydroxy-2-oxopropyl)-2-hexadecenoate REASON: '
               'MISSED No proper alcohol group found\n'
               ' * SMILES: '
               'O=C1O[C@@H](CCCCCCCCCCCCC[C@H](O)C)[C@H](C1=C)C(=O)O NAME: '
               'Allo-murolic acid REASON: MISSED No proper alcohol group '
               'found\n'
               ' * SMILES: OC(CCCCCCCCCCC)CC(O)CCC1=CC(OC)=C(OC)C=C1 NAME: '
               'Methyl-[12]-gingerdiol REASON: MISSED No proper alcohol group '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No proper alcohol group found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No proper alcohol group found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No proper alcohol group found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Not enough total carbons for a '
                                               'long hydrocarbon chain'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'Not enough total carbons for a '
                                               'long hydrocarbon chain'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'Not enough total carbons for a '
                                               'long hydrocarbon chain'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=C(C=C3)Cl)C(=O)N(C[C@@H]1OC)C)C)C',
                                     'name': '1-(4-chlorophenyl)-3-[(4S,7R,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'No proper alcohol group found'},
                                 {   'smiles': '[H][C@@]1(CNc2nc(N)[nH]c(=O)c2N1)[C@@H](O)[C@H](C)O',
                                     'name': 'sapropterin',
                                     'reason': 'Not enough total carbons for a '
                                               'long hydrocarbon chain'},
                                 {   'smiles': 'CCCCN(CCCC)CCNC(=O)CNC(=O)C1=CC=CS1',
                                     'name': 'N-[2-[2-(dibutylamino)ethylamino]-2-oxoethyl]-2-thiophenecarboxamide',
                                     'reason': 'No proper alcohol group found'},
                                 {   'smiles': 'ClC1=C(OC2=C(OC)C(OC3=C(Cl)C=C(O)C=C3Cl)=CC(=C2)O)C(Cl)=CC(=C1)O',
                                     'name': 'Rubetauphelin B',
                                     'reason': 'No proper alcohol group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCCCCCC(C)O',
                                      'name': 'tridecan-2-ol',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'O[C@@H]([C@H](O)CC/C=C/CO)CC/C=C/C=C',
                                      'name': 'Harzianumol G',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'CCCCCCC(O)CCCCCC',
                                      'name': 'tridecan-7-ol',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'C\\C=C\\C#CC#CC#CC#C[C@@H](O)CO',
                                      'name': 'Dehydrosafynol',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'CCCCCCCCC(O)CCCC',
                                      'name': 'tridecan-5-ol',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'OC/C=C/C#CC#CC#CCCCC',
                                      'name': '(E)-2-Tridecene-4,6,8-triyn-1-ol',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'CC(C)CCCCCCCCCCO',
                                      'name': '11-methyldodecan-1-ol',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'O[C@H]([C@@H](O)CC/C=C/CO)CC/C=C/C=C',
                                      'name': 'Harzianumol H',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'CCCCCCCC(O)CCCCC',
                                      'name': 'tridecan-6-ol',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'},
                                  {   'smiles': 'CCCCCCCCCC(O)CCC',
                                      'name': 'tridecan-4-ol',
                                      'reason': 'Not enough total carbons for '
                                                'a long hydrocarbon chain'}],
    'attempt': 2,
    'reasoning': 'To understand what went wrong with the previous attempt, '
                 "let's analyze the code and the outcomes:\n"
                 '\n'
                 '1. **Alcohol Group Detection**: \n'
                 '   - The pattern `[CX4H2,CX4H3][OX2H]` checks for an oxygen '
                 'atom connected to a saturated carbon, which is appropriate '
                 'for detecting an alcohol group. However, it may be too '
                 'strict in excluding secondary and tertiary alcohols where '
                 'the carbon might not always have hydrogens.\n'
                 '\n'
                 '2. **Carbon Chain Detection**: \n'
                 '   - The pattern '
                 '`[CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0][CR0]` '
                 'is designed to match a long carbon chain but does not '
                 'account for unsaturations (like double or triple bonds) or '
                 'branching, which can still qualify as a long-chain fatty '
                 'alcohol.\n'
                 '\n'
                 '3. **False Negatives Analysis**:\n'
                 '   - Many molecules missed by the previous code have '
                 'unsaturations or substitutions in their carbon chain or have '
                 'the alcohol group not strictly at the end but within the '
                 'chain (including diols). This suggests a need to broaden the '
                 'alcohol detection and allow for more flexible chain '
                 'structures.\n'
                 '\n'
                 '### Improvements:\n'
                 '\n'
                 '1. **Broaden Alcohol Pattern**: \n'
                 '   - Allow for primary, secondary, and tertiary alcohols by '
                 'adjusting the alcohol detection pattern.\n'
                 '\n'
                 '2. **Flexible Carbon Chain**: \n'
                 '   - Enable matching for longer chains that include '
                 'unsaturations and possibly branching. Ensure there are at '
                 'least 14 total carbons.\n'
                 '\n'
                 '3. **Combinatory Detection**:\n'
                 '   - Use comprehensive substructure searches and combined '
                 'rules to confirm both a long carbon chain and the presence '
                 'of an alcohol group.\n'
                 '\n'
                 "Here's how you might adapt the function:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 294,
    'num_false_positives': 100,
    'num_true_negatives': 98,
    'num_false_negatives': 12,
    'num_negatives': None,
    'precision': 0.7461928934010152,
    'recall': 0.9607843137254902,
    'f1': 0.8400000000000001,
    'accuracy': 0.7777777777777778,
    'negative_predictive_value': 0.8909090909090909}