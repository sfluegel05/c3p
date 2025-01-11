"""
Classifies: CHEBI:16180 N-acylglycine
"""
from rdkit import Chem

def is_N_acylglycine(smiles: str):
    """
    Determines if a molecule is an N-acylglycine based on its SMILES string.
    An N-acylglycine consists of a glycine moiety with an acyl group bonded to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an N-acylglycine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycine core pattern (more flexible)
    glycine_pattern = Chem.MolFromSmarts("NCC(=O)[O,N]")  # allows flexibility with acids or derivatives
    if not mol.HasSubstructMatch(glycine_pattern):
        return False, "No glycine core structure found"
    
    # N-acyl linkage pattern (specific to an acyl group attached via nitrogen)
    acyl_glycine_pattern = Chem.MolFromSmarts("C(=O)NCC(=O)[O,N]")  # carbonyl group bonded to nitrogen with glycine motif
    if not mol.HasSubstructMatch(acyl_glycine_pattern):
        return False, "No acyl group bonded to glycine nitrogen as expected for N-acylglycine"

    return True, "Contains N-acylglycine characteristic structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16180',
                          'name': 'N-acylglycine',
                          'definition': 'An N-acyl-amino acid in which amino '
                                        'acid specified is glycine.',
                          'parents': [   'CHEBI:140325',
                                         'CHEBI:24373',
                                         'CHEBI:51569'],
                          'xrefs': ['KEGG:C02055', 'MetaCyc:CPD-426'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: C(CNC(=O)CCC)(=O)O NAME: '
               'butyrylglycine REASON: MISSED No acyl group bonded to glycine '
               'nitrogen as expected for N-acylglycine\n'
               ' * SMILES: Cc1cccc(c1)C(=O)NCC(O)=O NAME: m-methylhippuric '
               'acid REASON: MISSED No acyl group bonded to glycine nitrogen '
               'as expected for N-acylglycine\n'
               ' * SMILES: '
               'ClC1=C(O)C(Cl)=CC(=C1)C(=O)NC(C(=O)O)CC2NC(C(=O)O)(C(O)C(N)C(=O)O)CC2 '
               'NAME: Kaitocephalin REASON: MISSED No acyl group bonded to '
               'glycine nitrogen as expected for N-acylglycine\n'
               ' * SMILES: FC/C(=N\\CCC[C@H](NC(=O)C1=CC=CC=C1)C(=O)N)/N NAME: '
               'F-Amidine REASON: MISSED No glycine core structure found\n'
               ' * SMILES: ClC=1C=C(C(=O)NCC(O)=O)C=CC1 NAME: '
               'm-Chloro-hippuric acid REASON: MISSED No acyl group bonded to '
               'glycine nitrogen as expected for N-acylglycine\n'
               ' * SMILES: O=C(NCC(=O)N/C=C/C=1C2=C(C=CC=C2)NC1)C3=C(N)C=CC=C3 '
               'NAME: Penidiamide REASON: MISSED No glycine core structure '
               'found\n'
               ' * SMILES: O=C(N)NCC(=O)O NAME: N-carbamoylglycine REASON: '
               'MISSED No acyl group bonded to glycine nitrogen as expected '
               'for N-acylglycine\n'
               ' * SMILES: OC(=O)CNC(=O)c1cccnc1 NAME: N-nicotinoylglycine '
               'REASON: MISSED No acyl group bonded to glycine nitrogen as '
               'expected for N-acylglycine\n'
               ' * SMILES: FC(F)(F)C=1C=C(C=CC1)C(=O)NCC(O)=O NAME: '
               'm-Trifluoromethylhippuric acid REASON: MISSED No acyl group '
               'bonded to glycine nitrogen as expected for N-acylglycine\n'
               ' * SMILES: OC(=O)CNC(=O)CCCP(O)(=O)Oc1ccc(cc1)[N+]([O-])=O '
               'NAME: N-(p-nitrophenylphosphobutanoyl)glycine REASON: MISSED '
               'No acyl group bonded to glycine nitrogen as expected for '
               'N-acylglycine\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@@H](C(=O)[C@@]2([C@]1(C(=C)C[C@H]3[C@]4([C@@H](C(C(=O)CC4)(C)C)[C@H](C[C@]23C)O)C)C)C(=O)OC)C',
                                     'name': 'Terreustoxin L',
                                     'reason': 'No glycine core structure '
                                               'found'},
                                 {   'smiles': '[H][C@]12CCC[C@H](C)[C@@]1(C)C[C@@H](CC2)C(C)C',
                                     'name': 'eremophilane',
                                     'reason': 'No glycine core structure '
                                               'found'},
                                 {   'smiles': 'N.C[C@]12CC[C@@](C)(C[C@H]1C1=CC(=O)[C@@H]3[C@@]4(C)CC[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)C(O)=O)C(O)=O)C(C)(C)[C@@H]4CC[C@@]3(C)[C@]1(C)CC2)C(O)=O',
                                     'name': 'Monoammonium glycyrrhizinate',
                                     'reason': 'No glycine core structure '
                                               'found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@H]2[C@H](O[C@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O)O[C@@H]1CO)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)[C@@H]4O[C@@H]6OC[C@@H](O)[C@H](O)[C@H]6O)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[(2R,3S,4R,5R,6R)-5-Acetamido-6-hydroxy-2-(hydroxymethyl)-4-[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-4-hydroxy-5-[(2S,3S,4S,5R,6R)-5-hydroxy-6-(hydroxymethyl)-4-[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-3-[(2S,3R,4S,5R)-3,4,5-trihydroxyoxan-2-yl]oxyoxan-2-yl]oxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No glycine core structure '
                                               'found'},
                                 {   'smiles': 'CN(C)CC(=O)N[C@H]1CC[C@@H](O[C@H]1CO)CCN2C=C(N=N2)C3=CC(=CC=C3)OC',
                                     'name': '2-(dimethylamino)-N-[(2R,3S,6R)-2-(hydroxymethyl)-6-[2-[4-(3-methoxyphenyl)-1-triazolyl]ethyl]-3-oxanyl]acetamide',
                                     'reason': 'No acyl group bonded to '
                                               'glycine nitrogen as expected '
                                               'for N-acylglycine'},
                                 {   'smiles': 'CC1=CC=C(C=C1)CS(=O)(=O)C2=NC=C(C(=N2)C(=O)NC3=NN=C(S3)C(C)C)Cl',
                                     'name': '5-chloro-2-[(4-methylphenyl)methylsulfonyl]-N-(5-propan-2-yl-1,3,4-thiadiazol-2-yl)-4-pyrimidinecarboxamide',
                                     'reason': 'No glycine core structure '
                                               'found'},
                                 {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\CCCCC)([O-])=O',
                                     'name': 'PC(15:1(9Z)/15:1(9Z))',
                                     'reason': 'No glycine core structure '
                                               'found'},
                                 {   'smiles': 'O=C(OC[C@]1(O)[C@@H]2C=C3C=C[C@@H]4[C@@]([C@]3(C2)CC1)(CC[C@H]([C@]4(CO)C)O)C)C',
                                     'name': 'Aphidicolin A42',
                                     'reason': 'No glycine core structure '
                                               'found'},
                                 {   'smiles': 'O=C1OC(=CC(=C1)O)[C@H]2C3=C(C(O)=CC=C3)C(=O)C[C@]2(O)C4=CC=CC=C4',
                                     'name': 'Wailupemycin D',
                                     'reason': 'No glycine core structure '
                                               'found'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C3=CC=CC=C3O2)NC(=O)CSC4=NC5=CC=CC=C5N4',
                                     'name': '2-(1H-benzimidazol-2-ylthio)-N-(2-methoxy-3-dibenzofuranyl)acetamide',
                                     'reason': 'No glycine core structure '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Upon analyzing the previous code and the outcomes, it seems '
                 'the logic for identifying the N-acylglycine structure did '
                 'not sufficiently handle variations in molecule structures, '
                 'especially considering complex branched patterns or diverse '
                 'functional groups. The main issues appear to be:\n'
                 '\n'
                 '1. The initial code for checking the existence of a glycine '
                 'core structure (NCC(=O)O) may not sufficiently capture all '
                 'instances due to the potential presence of substituents or '
                 'additional atoms.\n'
                 '2. The linkage detection between the glycine nitrogen and an '
                 'acyl group using the amide pattern (C(=O)N) might have '
                 'missed scenarios where atom indices are not placed where '
                 'expected.\n'
                 '\n'
                 'Improvements should include:\n'
                 '- More flexible pattern matching for the glycine component '
                 'to accommodate attached substituents.\n'
                 '- Ensuring the nitrogen in the glycine structure is '
                 'explicitly linked to an acyl (carbonyl group followed by a '
                 'chain of carbons) group.\n'
                 '\n'
                 "To increase the reliability of our pattern matching, here's "
                 'the revised approach:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 86,
    'num_false_positives': 100,
    'num_true_negatives': 724,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.46236559139784944,
    'recall': 1.0,
    'f1': 0.6323529411764706,
    'accuracy': 0.8901098901098901,
    'negative_predictive_value': 1.0}