"""
Classifies: CHEBI:73702 wax
"""
"""
Classifies: CHEBI:35195 wax
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_wax(smiles: str):
    """
    Determines if a molecule is a wax based on its SMILES string.
    A wax is typically a long-chain ester formed between a fatty acid and a fatty alcohol,
    characterized by long hydrocarbon chains and minimal other functional groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Must have exactly 1 ester group, found {len(ester_matches)}"

    # Count atoms and check composition
    total_atoms = mol.GetNumAtoms()
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    other_count = total_atoms - c_count - o_count
    
    # Waxes should be mostly carbon and hydrogen
    if c_count < 20:
        return False, f"Too few carbons ({c_count}) for a wax"
    if o_count != 2:
        return False, f"Waxes should have exactly 2 oxygens (ester group), found {o_count}"
    if other_count > 0:
        return False, "Contains atoms other than C, H, and O"

    # Check for cyclic structures (waxes are typically linear)
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0:
        return False, "Contains rings - waxes are typically linear"

    # Look for long carbon chains (at least 8 carbons in a row)
    long_chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    long_chain_matches = mol.GetSubstructMatches(long_chain_pattern)
    if len(long_chain_matches) < 2:
        return False, "Must have at least two long hydrocarbon chains"

    # Count rotatable bonds to verify flexibility
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Not enough rotatable bonds for a wax"

    # Check molecular weight - waxes typically >350 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 350:
        return False, "Molecular weight too low for a wax"
    if mol_wt > 1000:
        return False, "Molecular weight too high for typical wax"

    # Check for carboxylic acids
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if mol.HasSubstructMatch(acid_pattern):
        return False, "Contains free carboxylic acid group"

    # Calculate fraction of sp3 carbons (should be high for waxes)
    sp3_pattern = Chem.MolFromSmarts("[CX4]")
    sp3_count = len(mol.GetSubstructMatches(sp3_pattern))
    sp3_fraction = sp3_count / c_count
    if sp3_fraction < 0.8:
        return False, "Too many unsaturated carbons for a typical wax"

    return True, "Long-chain ester with appropriate molecular weight and composition typical of waxes"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73702',
                          'name': 'wax',
                          'definition': 'A chemical substance that is an '
                                        'organic compound or mixture of '
                                        'compounds that is composed of '
                                        'long-chain molecules and is malleable '
                                        'at ambient temperatures.',
                          'parents': ['CHEBI:59999', 'CHEBI:61697'],
                          'xrefs': [   'PMID:21333071',
                                       'PMID:23050609',
                                       'PMID:23103356',
                                       'PMID:23288373',
                                       'Wikipedia:Wax'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.01020408163265306 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CCCCCCCCCCCCCCCC(=O)OCC(CCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC '
               'NAME: 2-palmitoyloxypalmityl palmitate REASON: CORRECT '
               'Long-chain ester with appropriate molecular weight and chain '
               'length\n'
               'False positives: SMILES: '
               'O[C@H](CCCCCCCCCCCCC1=C[C@@H](OC1=O)C)[C@H](O)CC/C=C\\CCCCCCCCCCCC '
               'NAME: Cohibin A REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'C\\C=C(\\C)C(=O)O[C@H](C\\C(C)=C\\CO)[C@@H]1[C@H](OC(=O)C1=C)\\C=C(/C)CO '
               'NAME: Pycnolide REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: CC1C=C(CCCCCCCCCCCCCCCCC(C(CCCCCCCCCCCC)O)O)C(=O)O1 '
               'NAME: artemoin B REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               '[H][C@@]12CC(=O)O[C@]1(C)C[C@](C)(CCCCCCCCCCCCc1ccc(O)cc1)OO2 '
               'NAME: plakortolide M REASON: WRONGLY CLASSIFIED Long-chain '
               'ester with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'O=C1C=C[C@H]([C@@]1(O)C(OC(=O)C)CCCCCCCCCCCC)OC(=O)C NAME: '
               '4,6-di-O-acetyl hygrophorone A12 REASON: WRONGLY CLASSIFIED '
               'Long-chain ester with appropriate molecular weight and chain '
               'length\n'
               ' * SMILES: O=C1C=C[C@H]([C@]1(O)C(=O)CCCCCCCCCCCC)OC(=O)C '
               'NAME: 4-O-acetyl hygrophorone C12 REASON: WRONGLY CLASSIFIED '
               'Long-chain ester with appropriate molecular weight and chain '
               'length\n'
               ' * SMILES: '
               'O=C1C=C[C@@H]([C@@]1(OC(=O)C)C(OC(=O)C)CCCCCCCCCCCCCC)OC(=O)C '
               'NAME: 4,5,6-Tri-O-acetyl hygrophorone B14 REASON: WRONGLY '
               'CLASSIFIED Long-chain ester with appropriate molecular weight '
               'and chain length\n'
               ' * SMILES: '
               'C/C(=C\\CC\\C(\\C)=C\\C[C@]12C(C=C([C@H]([C@H]1O2)OC(=O)C)OC)=O)/CCC=C(C)C '
               'NAME: yanuthone X1 REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'C(CCCCCCCC(=O)OC(CCCCCCCCCCCC(=O)[O-])CCCCC)CCCCCCCCC NAME: '
               '13-(octadecanoyloxy)octadecanoate REASON: WRONGLY CLASSIFIED '
               'Long-chain ester with appropriate molecular weight and chain '
               'length\n'
               ' * SMILES: '
               'BrC1=C2C(C(=O)[C@](C1=O)(OC(=O)/C(=C/[C@H](CCCCCC)C)/C)C)=COC(=C2)C3=C(O)C=CC=C3C '
               'NAME: Penicilone H REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'O1C=2C(CC=C(C)C)(CC=C(C)C)C(=O)C(=C(O)C2C(CC1=O)C3=CC=CC=C3)C(=O)C(CC)C '
               'NAME: Calofloride REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'O=C1C=C2C(=COC(=C2)C)C=3[C@]1(OC(=O)C3C(=O)[C@H](CCCCCCCC)C)C '
               'NAME: Deflectin 2a REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'O=C1N(CC(=O)O[C@H]([C@@H](C(=O)[C@@](O)(C)C(NCC(NC(C=C1)=C)=O)=O)C)CCCCCCCCCCCCCCC)C '
               'NAME: BE-43547-C2 REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'O=C1O[C@](OCC)(C2=C(O)C=CC(=C2)O)C=C1CC/C=C(/CC/C=C(\\CO)/C)\\C '
               'NAME: Australeol E REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: O=C(C1(O)C(OC(=O)C)C=C[C@H]1OC(=O)C)CCCCCCCCCC '
               'NAME: 1,4-di-O-acetyl hygrophorone E10 REASON: WRONGLY '
               'CLASSIFIED Long-chain ester with appropriate molecular weight '
               'and chain length\n'
               ' * SMILES: '
               'O=C(OC(C(O)CCCCCCCC)CCCCCCCCCCCCCCCC1=CC(O)=CC(=C1)O)C NAME: '
               'SNF-4794-7 REASON: WRONGLY CLASSIFIED Long-chain ester with '
               'appropriate molecular weight and chain length\n'
               ' * SMILES: CCCCCCCCCCCC(CC1=CC(=O)C=C(OC)C1=O)OC(C)=O NAME: '
               'Cornudentanone REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'ClC1=C2C(C3=C(C(=O)CC)C(O[C@@]3(C1=O)C)=O)=CN(C/C=C(/CCC=C(C)C)\\C)C(=C2)/C=C/[C@H](CC)C '
               'NAME: '
               '4′,5′-dinor-5′-deoxy-N-(3,7-dimethyl-2,6-octadienyl)-2-aza-2-deoxychaetoviridin '
               'A REASON: WRONGLY CLASSIFIED Long-chain ester with appropriate '
               'molecular weight and chain length\n'
               ' * SMILES: '
               'O1C(C1CC/C=C\\CCCCCCCCCCCCC2=CC(OC2=O)C)CCCCCCCCCCCC NAME: '
               'Epoxymurin B REASON: WRONGLY CLASSIFIED Long-chain ester with '
               'appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'O=C1OC(C2=C(O)C=CC(=C2)O)C=C1CC/C=C(/C(O)C/C=C(/COC(=O)/C=C/C3=CC=C(O)C=C3)\\C)\\CO '
               'NAME: Zizhine H REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'O=C1C(=CCCC(C)=C[C@H](C(C=C1)(C)C)OC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)C '
               'NAME: Mitibetaimol A linoleate REASON: WRONGLY CLASSIFIED '
               'Long-chain ester with appropriate molecular weight and chain '
               'length\n'
               ' * SMILES: '
               'CC(=O)O[C@@H](C\\C(C)=C\\CCc1ccoc1)\\C=C(/C)CC[C@@H](OC(C)=O)[C@]1(C)CCC(=O)O1 '
               'NAME: irciformonin F REASON: WRONGLY CLASSIFIED Long-chain '
               'ester with appropriate molecular weight and chain length\n'
               ' * SMILES: CCCCCCCCCCCCCCCC(=O)OC(CCCCCCCCCCCCC)CCCC([O-])=O '
               'NAME: 5-PAHSA(1-) REASON: WRONGLY CLASSIFIED Long-chain ester '
               'with appropriate molecular weight and chain length\n'
               ' * SMILES: '
               'C(\\CCCCCCCC(=O)OC(CCCCCCCC(=O)[O-])CCCCCCCCC)=C\\CCCCCCCC '
               'NAME: 9-[(9Z)-octadecenoyloxy]octadecanoate REASON: WRONGLY '
               'CLASSIFIED Long-chain ester with appropriate molecular weight '
               'and chain length\n'
               ' * SMILES: '
               'C[C@](O)(CCCCCCCC\\C=C\\C=C\\c1ccccc1)C[C@@]1(C)OC(=O)C[C@H]1O '
               'NAME: seco-plakortolide O REASON: WRONGLY CLASSIFIED '
               'Long-chain ester with appropriate molecular weight and chain '
               'length\n'
               'False negatives: SMILES: '
               'O(CCCCCCCC/C=C\\CCCCCCCC)C(=O)CCCCCCCCCCCCCCCCCCCCC NAME: '
               'Oleyl behenate REASON: MISSED Simple methyl/ethyl ester, not a '
               'wax\n'
               ' * SMILES: O(CCCCCCCC/C=C\\CCCCCC)C(=O)CCCCCCCCCCCCCCCCCCCCC '
               'NAME: Palmitoleyl behenate REASON: MISSED Simple methyl/ethyl '
               'ester, not a wax\n'
               ' * SMILES: '
               'O(CCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(=O)CCCCCCC/C=C\\CCCCCCCC '
               'NAME: Linolenyl oleate REASON: MISSED Simple methyl/ethyl '
               'ester, not a wax\n'
               ' * SMILES: '
               'O(CCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(=O)CCCCCCC/C=C\\CCCCCC '
               'NAME: Linolenyl palmitatoleate REASON: MISSED Simple '
               'methyl/ethyl ester, not a wax\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCCCCOC(=O)CCCCCCC\\C=C/CCCCCCCC '
               'NAME: arachidyl oleate REASON: MISSED Simple methyl/ethyl '
               'ester, not a wax\n'
               ' * SMILES: O(CCCCCCCC/C=C\\CCCCCC)C(=O)CCCCCCC/C=C\\CCCC NAME: '
               'Palmitoleyl myristoleate REASON: MISSED Simple methyl/ethyl '
               'ester, not a wax\n'
               ' * SMILES: '
               'O(CCCCCCCC/C=C\\CCCCCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC '
               'NAME: Oleyl linolenate REASON: MISSED Simple methyl/ethyl '
               'ester, not a wax\n'
               ' * SMILES: '
               'O(CCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC '
               'NAME: Linoleyl linoleate REASON: MISSED Simple methyl/ethyl '
               'ester, not a wax\n'
               ' * SMILES: O(CCCCCCCC/C=C\\CCCCCC)C(=O)CCCCCCCCCCCCCCC NAME: '
               'Palmitoleyl palmitate REASON: MISSED Simple methyl/ethyl '
               'ester, not a wax\n'
               ' * SMILES: O(CCCCCCCCCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC NAME: '
               'Lauryl linoleate REASON: MISSED Simple methyl/ethyl ester, not '
               'a wax\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC(=O)CCNC(O)=O',
                                     'name': 'N-carboxy-beta-alanine',
                                     'reason': 'Must have exactly 1 ester '
                                               'group, found 2'},
                                 {   'smiles': 'C[C@@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)NC3=CC=CC=C3F)C(=O)N(C[C@H]1OC)C)C)CCC(F)(F)F',
                                     'name': '1-(2-fluorophenyl)-3-[(5S,6R,9R)-5-methoxy-3,6,9-trimethyl-2-oxo-8-(3,3,3-trifluoropropyl)-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'Must have exactly 1 ester '
                                               'group, found 0'},
                                 {   'smiles': 'O(C1=CC(=C(O)C=C1)C2=C(O)C=CC(=C2)OC)C',
                                     'name': '5,5′-dimethoxybiphenyl-2,2′-diol',
                                     'reason': 'Must have exactly 1 ester '
                                               'group, found 0'},
                                 {   'smiles': 'COC1=CC=C(CCNC[C@H](O)C2=CC=C(O)C=C2)C=C1OC',
                                     'name': 'denopamine',
                                     'reason': 'Must have exactly 1 ester '
                                               'group, found 0'},
                                 {   'smiles': '[Zn++].[S-]C(=S)NCCNC([S-])=S',
                                     'name': 'zineb',
                                     'reason': 'Must have exactly 1 ester '
                                               'group, found 0'},
                                 {   'smiles': 'P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCCCCCCC(C)C)COC(=O)CCCCCCCCC(CC)C)(O)=O)(O)=O',
                                     'name': 'CDP-DG(a-13:0/i-15:0)',
                                     'reason': 'Must have exactly 1 ester '
                                               'group, found 2'},
                                 {   'smiles': 'O=C1OC(O)C2=C1C[C@](O)([C@H]3CC(C[C@H]3[C@@H]2O)(C)C)C',
                                     'name': 'Lactarolide A',
                                     'reason': 'Too few carbons (15) for a '
                                               'wax'},
                                 {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'Too few carbons (16) for a '
                                               'wax'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO',
                                     'name': 'N-[(3R,4R,5S,6R)-2,5-Dihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-3-yl]acetamide',
                                     'reason': 'Must have exactly 1 ester '
                                               'group, found 0'},
                                 {   'smiles': 'COc1ccc(\\C=C/C2CCC=CC2c2ccc(OC)c(OC)c2)cc1OC',
                                     'name': '3-(3,4-dimethoxyphenyl)-4-[(Z)-3,4-dimethoxystyryl]cyclohex-1-ene',
                                     'reason': 'Must have exactly 1 ester '
                                               'group, found 0'}],
    'sample_false_negatives': [   {   'smiles': 'O(CCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(=O)CCCCCCC/C=C\\CCCCCCCC',
                                      'name': 'Linolenyl oleate',
                                      'reason': 'Too many unsaturated carbons '
                                                'for a typical wax'},
                                  {   'smiles': 'O(CCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(=O)CCCCCCC/C=C\\CCCCCC',
                                      'name': 'Linolenyl palmitatoleate',
                                      'reason': 'Too many unsaturated carbons '
                                                'for a typical wax'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OCC(CCCCCCCCCCCCCC)OC(=O)CCCCCCCCCCCCCCC',
                                      'name': '2-palmitoyloxypalmityl '
                                              'palmitate',
                                      'reason': 'Must have exactly 1 ester '
                                                'group, found 2'},
                                  {   'smiles': 'O(CCCCCCCC/C=C\\CCCCCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC',
                                      'name': 'Oleyl linolenate',
                                      'reason': 'Too many unsaturated carbons '
                                                'for a typical wax'},
                                  {   'smiles': 'O(CCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC',
                                      'name': 'Linoleyl linoleate',
                                      'reason': 'Too many unsaturated carbons '
                                                'for a typical wax'},
                                  {   'smiles': 'O(CCCCCCCC/C=C\\C/C=C\\CCCCC)C(=O)CCCCCCC/C=C\\CCCCCC',
                                      'name': 'Linoleyl palmitatoleate',
                                      'reason': 'Too many unsaturated carbons '
                                                'for a typical wax'},
                                  {   'smiles': 'O(CCCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(=O)CCCCCCCCCCCCC',
                                      'name': 'Linolenyl myristate',
                                      'reason': 'Too many unsaturated carbons '
                                                'for a typical wax'},
                                  {   'smiles': 'OC(CCCCCCCCOC(=O)CCCCCCCCCCCCCCCCCCCCC)CCCC',
                                      'name': '9-Hydroxytridecyl docosanoate',
                                      'reason': 'Waxes should have exactly 2 '
                                                'oxygens (ester group), found '
                                                '3'},
                                  {   'smiles': 'O(CCCCCCCC/C=C\\CCCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\CCCCC',
                                      'name': 'Palmitoleyl linoleate',
                                      'reason': 'Too many unsaturated carbons '
                                                'for a typical wax'},
                                  {   'smiles': 'O(CCCCCCCC/C=C\\CCCCCC)C(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC',
                                      'name': 'Palmitoleyl linolenate',
                                      'reason': 'Too many unsaturated carbons '
                                                'for a typical wax'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous outcomes and suggest '
                 'improvements:\n'
                 '\n'
                 '1. MAIN ISSUE: The program failed to identify many '
                 'legitimate waxes (false negatives) while incorrectly '
                 'classifying many non-waxes as waxes (false positives).\n'
                 '\n'
                 '2. ANALYSIS OF FALSE NEGATIVES:\n'
                 '- All missed compounds were legitimate wax esters with long '
                 'hydrocarbon chains\n'
                 '- The pattern shows these are mostly fatty acid esters with '
                 'fatty alcohols\n'
                 '- The check for "ester group not connected to long chains on '
                 'both sides" might be too strict\n'
                 '\n'
                 '3. ANALYSIS OF FALSE POSITIVES:\n'
                 '- Many complex natural products with ester groups were '
                 'incorrectly classified\n'
                 '- Many cyclic compounds were mistakenly identified as waxes\n'
                 '- Compounds with multiple functional groups beyond simple '
                 'esters were misclassified\n'
                 '\n'
                 '4. IMPROVEMENTS NEEDED:\n'
                 '- Remove overly restrictive chain connectivity checks\n'
                 '- Add checks to exclude cyclic structures (waxes are '
                 'typically linear)\n'
                 '- Be more specific about the required chain lengths\n'
                 '- Add checks to exclude compounds with too many '
                 'non-carbon/hydrogen atoms\n'
                 '- Focus on the ratio of carbons to other atoms (waxes are '
                 'mostly hydrocarbons)\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 73,
    'num_false_positives': 21,
    'num_true_negatives': 142184,
    'num_false_negatives': 22,
    'num_negatives': None,
    'precision': 0.776595744680851,
    'recall': 0.7684210526315789,
    'f1': 0.7724867724867724,
    'accuracy': 0.9996978215038651,
    'negative_predictive_value': 0.9998452948539442}