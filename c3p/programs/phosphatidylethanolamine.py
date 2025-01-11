"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is a glycerophospholipid where a phosphatidyl group is esterified to
    the hydroxy group of ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Glycerol backbone with correct stereochemistry (chiral center at C2)
    glycerol_pattern = Chem.MolFromSmarts("[C@H](CO[P](=O)(O)OCCN)OC(=O)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone with correct stereochemistry"

    # Look for a second ester linkage on the other carbon of glycerol
    second_ester_pattern = Chem.MolFromSmarts("OC(=O)C")
    ester_matches = mol.GetSubstructMatches(second_ester_pattern)
    if len(ester_matches) < 2:
        return False, "Missing second ester linkage, found less than 2"

    return True, "Structure matches phosphatidylethanolamine class"

# Example usage
smiles_example = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OCCN)(O)=O"
result, reason = is_phosphatidylethanolamine(smiles_example)
print(result, reason)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16038',
                          'name': 'phosphatidylethanolamine',
                          'definition': 'A class of glycerophospholipids in '
                                        'which a phosphatidyl group is '
                                        'esterified to the hydroxy group of '
                                        'ethanolamine.',
                          'parents': ['CHEBI:36314'],
                          'xrefs': [   'DrugBank:DB04327',
                                       'HMDB:HMDB0060501',
                                       'KEGG:C00350',
                                       'LIPID_MAPS_instance:LMGP02010000',
                                       'PMID:10540156',
                                       'PMID:11042504',
                                       'PMID:11159918',
                                       'PMID:11829744',
                                       'PMID:12139474',
                                       'PMID:15653902',
                                       'PMID:16037249',
                                       'PMID:16303767',
                                       'PMID:16620109',
                                       'PMID:18034796',
                                       'PMID:18259190',
                                       'PMID:18398168',
                                       'PMID:18462396',
                                       'PMID:18570887',
                                       'PMID:18957134',
                                       'PMID:19393163',
                                       'PMID:23354482',
                                       'PMID:23369752',
                                       'PMID:23543734',
                                       'PMID:3196084',
                                       'PMID:7980848',
                                       'Wikipedia:Phosphatidylethanolamine'],
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
               'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCC)(OCCN)(O)=O '
               'NAME: PE(18:1(11Z)/16:0) REASON: MISSED No glycerol backbone '
               'found with correct stereochemistry\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCC)(OCCN)(O)=O '
               'NAME: PE(14:1(9Z)/17:0) REASON: MISSED No glycerol backbone '
               'found with correct stereochemistry\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O '
               'NAME: PE(18:1(9Z)/22:5(7Z,10Z,13Z,16Z,19Z)) REASON: MISSED No '
               'glycerol backbone found with correct stereochemistry\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(OCCN)(O)=O '
               'NAME: PE(18:2(9Z,12Z)/17:0) REASON: MISSED No glycerol '
               'backbone found with correct stereochemistry\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OCCN)(O)=O '
               'NAME: PE(17:1(9Z)/17:1(9Z)) REASON: MISSED No glycerol '
               'backbone found with correct stereochemistry\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCC(O)=O)COC(=O)CCCCCCCCCCCCCCC)(OCCN)(O)=O '
               'NAME: PG-PE REASON: MISSED No glycerol backbone found with '
               'correct stereochemistry\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O '
               'NAME: PE(20:1(11Z)/20:5(5Z,8Z,11Z,14Z,17Z)) REASON: MISSED No '
               'glycerol backbone found with correct stereochemistry\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OCCN)OC(=O)CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC '
               'NAME: '
               '1-octadecanoyl-2-(7Z,10Z,13Z,16Z)-docosatetraenoyl-sn-glycero-3-phosphoethanolamine '
               'REASON: MISSED No glycerol backbone found with correct '
               'stereochemistry\n'
               ' * SMILES: '
               'P(OCC(OC(=O)CCCCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\CCCCC)(OCCNC)(O)=O '
               'NAME: PE-NMe(18:2(9Z,12Z)/18:0) REASON: MISSED No glycerol '
               'backbone found with correct stereochemistry\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O '
               'NAME: PE(18:3(9Z,12Z,15Z)/20:2(11Z,14Z)) REASON: MISSED No '
               'glycerol backbone found with correct stereochemistry\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No glycerol backbone with '
                                               'correct stereochemistry'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCC\\C=C/CCCCCCCC(=O)OCC(O)COP(O)(=O)OCCN',
                                      'name': '1-oleoyl '
                                              'phosphatidylethanolamine',
                                      'reason': 'No glycerol backbone with '
                                                'correct stereochemistry'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](O)COC(=O)CCCCCCCC)([O-])=O',
                                      'name': 'PC(9:0/0:0)',
                                      'reason': 'No glycerol backbone with '
                                                'correct stereochemistry'}],
    'attempt': 1,
    'reasoning': 'To classify phosphatidylethanolamine (PE) molecules '
                 'correctly, we need to better analyze the structure of PE. A '
                 'phosphatidylethanolamine generally consists of:\n'
                 '\n'
                 '1. **Glycerol Backbone**: The central core is a glycerol '
                 'molecule. The correct stereochemistry at the carbon centers '
                 'should be captured.\n'
                 '\n'
                 '2. **Phosphate Group**: A phosphate group attached to the '
                 'glycerol backbone.\n'
                 '\n'
                 '3. **Ethanolamine Group**: An ethanolamine group is '
                 'esterified to the phosphate.\n'
                 '\n'
                 '4. **Two Fatty Acid Chains**: Typically attached to the '
                 'glycerol backbone through ester linkage, represented as '
                 '`[O-C(=O)-C...C]`.\n'
                 '\n'
                 'Given the issues with the previous code, the following '
                 'modifications are necessary:\n'
                 '\n'
                 '- Ensure stereochemistry in the glycerol backbone is '
                 'accurately captured using SMARTS.\n'
                 '- The phosphate group detection should be integrated into '
                 "how it's attached to the glycerol.\n"
                 '- Confirm the presence of ethanolamine.\n'
                 '- Validate the presence of two fatty acid chains.\n'
                 '\n'
                 'With the above considerations, the following code defines '
                 'the `is_phosphatidylethanolamine` function:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 818,
    'num_false_positives': 100,
    'num_true_negatives': 12010,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.8910675381263616,
    'recall': 0.9975609756097561,
    'f1': 0.9413118527042577,
    'accuracy': 0.9921113689095128,
    'negative_predictive_value': 0.9998334998334998}