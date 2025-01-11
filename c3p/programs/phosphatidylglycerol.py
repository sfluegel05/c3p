"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
from rdkit import Chem

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern to detect the glycerol backbone part of a phosphatidylglycerol
    # Generalized for possible stereo/isomer configurations and common modifications
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO[P](=O)(O)O)(CO)O")  # Search for glycerol-phosphate linkage
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)

    if not glycerol_matches:
        return False, "No glycerol backbone with phosphate linkage found"

    # Pattern for ester-linked long fatty acid chains
    ester_fatty_acid_pattern = Chem.MolFromSmarts("OC(=O)C")  # This captures ester linkage to a fatty acid chain
    ester_matches = mol.GetSubstructMatches(ester_fatty_acid_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester-linked fatty acids, found {len(ester_matches)}"

    # Ensure presence of at least one phosphatidyl group
    phosphatidyl_pattern = Chem.MolFromSmarts("COP(O)(=O)OC")
    phosphatidyl_matches = mol.GetSubstructMatches(phosphatidyl_pattern)
    if not phosphatidyl_matches:
        return False, "Missing phosphatidyl group"

    return True, "Contains glycerol backbone with phosphatidyl group and ester-linked fatty acid chains"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17517',
                          'name': 'phosphatidylglycerol',
                          'definition': 'A glycerophosphoglycerol that is '
                                        'glycerol in which the hydrogen of one '
                                        'of the primary hydroxy groups has '
                                        'been replaced by a phosphatidyl '
                                        'group.',
                          'parents': ['CHEBI:24360'],
                          'xrefs': ['KEGG:C00344'],
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
               'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O '
               'NAME: PG(18:1(11Z)/18:1(11Z)) REASON: MISSED No '
               'phosphatidylglycerol backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCC(O)/C=C/C=O)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OC[C@@H](O)CO)(O)=O '
               'NAME: OHOOA-PG REASON: MISSED No phosphatidylglycerol backbone '
               'found\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCC(O)=O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O '
               'NAME: PS-PG REASON: MISSED No phosphatidylglycerol backbone '
               'found\n'
               ' * SMILES: '
               'C([C@@](COC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)(OC(CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O '
               'NAME: PG(18:2(9Z,12Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z)) REASON: '
               'MISSED No phosphatidylglycerol backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O '
               'NAME: PG(16:0/13:0) REASON: MISSED No phosphatidylglycerol '
               'backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O '
               'NAME: PG(12:0/20:5(5Z,8Z,11Z,14Z,17Z)) REASON: MISSED No '
               'phosphatidylglycerol backbone found\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@H](COP(O)(=O)OC[C@@H](O)CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC '
               "NAME: 1,2-dioleoyl-sn-glycero-3-phospho-(1'-sn-glycerol) "
               'REASON: MISSED No phosphatidylglycerol backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(OC[C@@H](O)CO)(O)=O '
               'NAME: PG(20:4(5Z,8Z,11Z,14Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z)) '
               'REASON: MISSED No phosphatidylglycerol backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O '
               'NAME: PG(22:0/18:4(6Z,9Z,12Z,15Z)) REASON: MISSED No '
               'phosphatidylglycerol backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCCCCCCC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)(OC[C@@H](O)CO)(O)=O '
               'NAME: PG(20:5(5Z,8Z,11Z,14Z,17Z)/13:0) REASON: MISSED No '
               'phosphatidylglycerol backbone found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'},
                                 {   'smiles': 'O1[C@H]([C@H](NC(=O)C)[C@@H](O)C[C@]1(O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O)C(O)=O)[C@H](O)[C@H](O[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](OC(=O)C)[C@H](O)CO)C(O)=O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-6-[(1S,2R)-2-[(2S,4S,5R,6R)-5-acetamido-6-[(1R,2R)-1-acetyloxy-2,3-dihydroxypropyl]-2-carboxy-4-hydroxyoxan-2-yl]oxy-1,3-dihydroxypropyl]-2-[(2R,3S,4S,5R,6S)-3,5-dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-4-hydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No glycerol backbone with '
                                               'phosphate linkage found'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Analyzing the previous attempt, it seems that the failure '
                 'might be due to the use of overly simplistic SMARTS patterns '
                 'which did not generalize well to the structures presented. '
                 'Specifically, the pattern for the glycerol backbone with a '
                 'phosphatidyl group might be too narrow and not matching the '
                 'diverse representation seen in the dataset.\n'
                 '\n'
                 'Key observations:\n'
                 '1. **Phosphatidyl Group:** The phosphatidylglycerol backbone '
                 'involves a glycerol unit connected via a phosphoester '
                 'linkage to another glycerol, with at least one being part of '
                 'a phosphatidyl structure.\n'
                 '2. **Ester Linkages:** These are essential for connecting '
                 'fatty acid chains. Two are necessary for defining this class '
                 'but their connection to glycerol needs refinement.\n'
                 '3. **Diverse Structures:** Some of the provided examples '
                 'have complex modifications such as additional sugars or '
                 "stereoisomers which weren't considered.\n"
                 '\n'
                 'To improve:\n'
                 '- Broaden the glycerol backbone detection pattern, focusing '
                 'not only on one hydroxyl group but providing greater '
                 'flexibility.\n'
                 '- Ensure correct identification of ester linkages connected '
                 'to a valid phosphatidylglycerol backbone. Both major ester '
                 'linkages (to fatty acids) must be included.\n'
                 '- Re-assess the extraction of phosphate groups. Verify that '
                 "the glycerol's hydroxyl is replaced appropriately.\n"
                 '\n'
                 'Here’s an updated version with these considerations:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 329,
    'num_false_positives': 100,
    'num_true_negatives': 5818,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.7668997668997669,
    'recall': 1.0,
    'f1': 0.8680738786279684,
    'accuracy': 0.9839923163118297,
    'negative_predictive_value': 1.0}