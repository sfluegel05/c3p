"""
Classifies: CHEBI:17088 monoacyl-sn-glycerol 3-phosphate
"""
"""
Classifies: CHEBI:57980 monoacyl-sn-glycerol 3-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_monoacyl_sn_glycerol_3_phosphate(smiles: str):
    """
    Determines if a molecule is a monoacyl-sn-glycerol 3-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoacyl-sn-glycerol 3-phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Create and verify SMARTS patterns
    patterns = {
        'phosphate': '[OX2,O-]P(=O)([OX2,O-])[OX2]',  # Phosphate group
        'ester': '[CX3](=O)[OX2]',  # Ester group
        'glycerol_backbone': '[CH2X4,CH2][CH,CHX4][CH2X4,CH2]',  # More flexible glycerol pattern
        'ethanolamine': 'OCCN',  # Pattern to exclude
        'choline': 'OCC[N+]',    # Pattern to exclude
        'ether': '[CH2,CH][OX2][CH2]C'  # Pattern to exclude ether lipids
    }
    
    compiled_patterns = {}
    for name, pattern in patterns.items():
        pat = Chem.MolFromSmarts(pattern)
        if pat is None:
            return False, f"Invalid SMARTS pattern for {name}"
        compiled_patterns[name] = pat

    # Exclude molecules with ethanolamine or choline groups
    if mol.HasSubstructMatch(compiled_patterns['ethanolamine']):
        return False, "Contains ethanolamine group"
    if mol.HasSubstructMatch(compiled_patterns['choline']):
        return False, "Contains choline group"

    # Check for phosphate group
    if not mol.HasSubstructMatch(compiled_patterns['phosphate']):
        return False, "No phosphate group found"

    # Check for single ester group (acyl chain)
    ester_matches = mol.GetSubstructMatches(compiled_patterns['ester'])
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for glycerol backbone
    if not mol.HasSubstructMatch(compiled_patterns['glycerol_backbone']):
        return False, "No glycerol backbone found"

    # Exclude ether lipids
    if mol.HasSubstructMatch(compiled_patterns['ether']):
        return False, "Contains ether linkage"

    # Count oxygen atoms - should be 7 
    # (3 from phosphate, 2 from ester, 2 from glycerol backbone)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count != 7:
        return False, f"Found {o_count} oxygen atoms, expected 7"

    # Count phosphorus atoms - should be exactly 1
    p_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
    if p_count != 1:
        return False, f"Found {p_count} phosphorus atoms, need exactly 1"

    # Additional check for acyl chain length (at least 4 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 7:  # minimum 3 for glycerol + 4 for shortest acyl chain
        return False, "Acyl chain too short"

    return True, "Contains glycerol backbone with one acyl chain and phosphate group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17088',
                          'name': 'monoacyl-sn-glycerol 3-phosphate',
                          'definition': 'An sn-glycero-3-phosphate compound '
                                        'having a single unspecified acyl '
                                        'group at either position 1 or '
                                        'position 2.',
                          'parents': ['CHEBI:16961', 'CHEBI:26706'],
                          'xrefs': ['KEGG:C03849'],
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
               'P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCC)(O)(O)=O NAME: '
               '1-(9Z-tetradecenoyl)-glycero-3-phosphate REASON: MISSED '
               'Incorrect number of oxygen atoms (7), expected 6\n'
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)(O)(O)=O NAME: '
               '1-(11Z-docosenoyl)-glycero-3-phosphate REASON: MISSED '
               'Incorrect number of oxygen atoms (7), expected 6\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)O[C@H](CO)COP(O)(O)=O '
               'NAME: 2-arachidonoyl-sn-glycero-3-phosphate REASON: MISSED '
               'Incorrect number of oxygen atoms (7), expected 6\n'
               ' * SMILES: P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCCCCC)(O)(O)=O '
               'NAME: PA(19:1(9Z)/0:0) REASON: MISSED Incorrect number of '
               'oxygen atoms (7), expected 6\n'
               ' * SMILES: P(OCC(O)COC(=O)CCCCCCCCCCCCCCCCCCCCC(CC)C)(O)(O)=O '
               'NAME: LysoPA(a-25:0/0:0) REASON: MISSED Incorrect number of '
               'oxygen atoms (7), expected 6\n'
               ' * SMILES: P(OCC(O)COC(=O)CCCCC)(O)(O)=O NAME: LPA 6:0 REASON: '
               'MISSED Incorrect number of oxygen atoms (7), expected 6\n'
               ' * SMILES: P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCCCC)(O)(O)=O '
               'NAME: 1-heneicosanoyl-glycero-3-phosphate REASON: MISSED '
               'Incorrect number of oxygen atoms (7), expected 6\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(O)=O '
               'NAME: 1-icosanoyl-sn-glycero-3-phosphate REASON: MISSED '
               'Incorrect number of oxygen atoms (7), expected 6\n'
               ' * SMILES: P(OCC(OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)CO)(O)(O)=O '
               'NAME: LysoPA(0:0/18:2(9Z,12Z)) REASON: MISSED Incorrect number '
               'of oxygen atoms (7), expected 6\n'
               ' * SMILES: P(OCC(OC(=O)CCCCCCC/C=C\\CCCCCCCC)CO)(O)(O)=O NAME: '
               'LysoPA(0:0/18:1(9Z)) REASON: MISSED Incorrect number of oxygen '
               'atoms (7), expected 6\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'S(=O)(=O)(C1=C2C(=CC3=C1NC=4C=CC=CC34)[C@@]5([C@H]([C@](C(=O)O)([C@@H](O)CC5)C)CC2)C)C6=CC7=C(NC8=C7C=C9[C@@]%10([C@H]([C@](C(=O)O)([C@@H](O)CC%10)C)CCC9=C8)C)C=C6',
                                     'name': 'Sulfadixiamycin C',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'CNC(O)=O',
                                     'name': 'methylcarbamic acid',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'CCNC(=O)NC1=CC2=C(C=C1)OC[C@H]3[C@@H](CC[C@H](O3)CC(=O)N[C@@H](C)C4=CC=CC=C4)N(C2=O)C',
                                     'name': '2-[(2S,4aR,12aR)-8-(ethylcarbamoylamino)-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'Contains ethanolamine group'},
                                 {   'smiles': 'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO',
                                     'name': 'CID 91851985',
                                     'reason': 'Contains ethanolamine group'},
                                 {   'smiles': 'O(C(=O)C(C1C(CN2C(C1)C=3NC=4C(C3CC2)=CC=CC4)CC)=COC)C',
                                     'name': 'Methyl '
                                             '2-(3-ethyl-1,2,3,4,6,7,12,12b-octahydroindolo[2,3-a]quinolizin-2-yl)-3-methoxyprop-2-enoate',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'O[C@H](/C=C/C=C/C=C/[C@H](O)[C@H](O)C=C)[C@H](O)/C=C/C',
                                     'name': 'Separacene C',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](O)C(O)O[C@H](CO)[C@@H]3O)[C@H](NC(C)=O)[C@H]2O[C@@H]2O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->4)-[beta-D-Galp3S-(1->3)]-beta-D-GlcpNAc-(1->3)-D-Galp',
                                     'reason': 'Contains ethanolamine group'},
                                 {   'smiles': 'C1=CC=CC2=C1C(N([C@H](C(N2)=O)CC=3C=CC(=CC3)OC)C)=O',
                                     'name': "(S)-4'-methoxycyclopeptine",
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C',
                                     'name': 'Val-Val-Trp',
                                     'reason': 'Contains ethanolamine group'},
                                 {   'smiles': 'C=1C(=C(C=CC1/C=C/CO)OC(CO)C(O)C=2C=C(C(=CC2)O)OC)OC',
                                     'name': 'guaiacylglycerol beta-coniferyl '
                                             'ether',
                                     'reason': 'No phosphate group found'}],
    'sample_false_negatives': [   {   'smiles': 'C(CCCCCCC/C=C\\CCCCCCCC)(=O)O[C@@H](COP(=O)(O)O)CO/C=C\\CCCCCCCCCCCCCCCC',
                                      'name': '1-[(1Z)-octadecenyl]-2-oleoyl-sn-glycero-3-phosphate',
                                      'reason': 'Contains ether linkage'},
                                  {   'smiles': 'CCCCCCCC\\C=C/CCCCCCCCOC[C@H](COP(O)(O)=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC',
                                      'name': '1-(9Z-octadecenyl)-2-arachidonoyl-sn-glycero-3-phosphate',
                                      'reason': 'Contains ether linkage'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCOC[C@H](COP(O)(O)=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC',
                                      'name': '1-palmityl-2-arachidonoyl-sn-glycero-3-phosphate',
                                      'reason': 'Contains ether linkage'}],
    'attempt': 3,
    'reasoning': 'Let me analyze the previous attempt and its failures:\n'
                 '\n'
                 '1. The key issue appears to be the oxygen count check. '
                 'Looking at the failed examples, they all contain 7 oxygens, '
                 'not 6 as we required. This makes sense because a '
                 'monoacyl-sn-glycerol 3-phosphate should have:\n'
                 '   - 3 oxygens in the phosphate group (P(=O)(OH)2)\n'
                 '   - 2 oxygens in the ester linkage (R-C(=O)-O-)\n'
                 '   - 2 oxygens in the glycerol backbone (one free OH and one '
                 'connecting to phosphate)\n'
                 '   Total: 7 oxygens, not 6 as previously specified\n'
                 '\n'
                 '2. The previous code was also too restrictive in other '
                 'ways:\n'
                 '   - The hydroxyl group check might be failing because some '
                 'hydroxyls could be deprotonated\n'
                 '   - The glycerol backbone pattern might be too specific\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 45,
    'num_false_positives': 24,
    'num_true_negatives': 142228,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.6521739130434783,
    'recall': 0.9375,
    'f1': 0.7692307692307693,
    'accuracy': 0.9998102600140548,
    'negative_predictive_value': 0.9999789075517995}