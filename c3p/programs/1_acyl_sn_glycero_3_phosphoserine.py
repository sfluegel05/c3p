"""
Classifies: CHEBI:52603 1-acyl-sn-glycero-3-phosphoserine
"""
from rdkit import Chem

def is_1_acyl_sn_glycero_3_phosphoserine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES string.
    References structures for this class typically include glycerophosphoserine backbone and
    an acyl substituent at the 1-hydroxy position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoserine, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the glycerophosphoserine backbone
    glycerophosphoserine_pattern = Chem.MolFromSmarts("P(OC[C@H](O)CO)(OC[C@H](N)C(=O)O)O")
    if not mol.HasSubstructMatch(glycerophosphoserine_pattern):
        return False, "No glycerophosphoserine backbone found"

    # Define the SMARTS pattern for an acyl group linked through an ester at the primary position
    acyl_hydroxy_pattern = Chem.MolFromSmarts("C(=O)OC[C@H](O)C")
    if not mol.HasSubstructMatch(acyl_hydroxy_pattern):
        return False, "No acyl group found at the 1-hydroxy position"
    
    # If both patterns are found, this is a positive match
    return True, "Contains 1-acyl-sn-glycero-3-phosphoserine structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52603',
                          'name': '1-acyl-sn-glycero-3-phosphoserine',
                          'definition': 'An sn-glycerophosphoserine compound '
                                        'having an acyl substituent at the '
                                        '1-hydroxy position.',
                          'parents': ['CHEBI:68510'],
                          'xrefs': ['PMID:16341241'],
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
               'P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(16:0/0:0) REASON: MISSED No glycerophosphoserine '
               'backbone found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OC[C@H](N)C(O)=O '
               'NAME: 1-stearoyl-sn-glycero-3-phosphoserine REASON: MISSED No '
               'glycerophosphoserine backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/0:0) REASON: MISSED No '
               'glycerophosphoserine backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(15:0/0:0) REASON: MISSED No glycerophosphoserine '
               'backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(16:1(9Z)/0:0) REASON: MISSED No glycerophosphoserine '
               'backbone found\n'
               ' * SMILES: '
               'P(OCC(O)COC(=O)CCCCCCC/C=C/CCCCCCCC)(OCC(N)C(O)=O)(O)=O NAME: '
               '1-Oleoylglycerophosphoserine REASON: MISSED No '
               'glycerophosphoserine backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCCCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(22:0/0:0) REASON: MISSED No glycerophosphoserine '
               'backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\CCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(14:1(9Z)/0:0) REASON: MISSED No glycerophosphoserine '
               'backbone found\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OC[C@H](N)C(O)=O '
               'NAME: 1-O-oleoyl-sn-glycero-3-phosphoserine REASON: MISSED No '
               'glycerophosphoserine backbone found\n'
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(22:1(11Z)/0:0) REASON: MISSED No glycerophosphoserine '
               'backbone found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No glycerophosphoserine '
                                               'backbone found'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'The previous code attempts to classify a chemical entity as '
                 'a 1-acyl-sn-glycero-3-phosphoserine based on its SMILES '
                 'string. It does so by looking for a specific '
                 'glycerophosphoserine backbone and an acyl group attached at '
                 'the primary hydroxy position. However, it failed to classify '
                 "the example structures correctly. Here's what might have "
                 'gone wrong and how we can improve the code:\n'
                 '\n'
                 '1. **Pattern Matching for Backbone**: The SMARTS pattern '
                 'used for detecting the glycerophosphoserine backbone may be '
                 'too specific or incorrect, missing structures that do have '
                 'this backbone. The pattern should be revised to accurately '
                 'encompass variations in orientation or slight structural '
                 'differences common in this class.\n'
                 '\n'
                 '2. **Matching Acyl Group at the Correct Position**: The way '
                 'the code checks for the acyl group at the 1-hydroxy position '
                 'might be flawed. It might not be correctly associating the '
                 'ester linkage to the right oxygen. We should ensure the '
                 'pattern unambiguously identifies the acyl linkage.\n'
                 '\n'
                 '3. **Stereo Centers**: Stereo configurations (indicated by @ '
                 'in SMILES strings) might be present and need consideration '
                 'for correct matching.\n'
                 '\n'
                 '4. **Validation Among Examples**: Ensure the created '
                 'patterns match all example structures to validate correct '
                 'detection, using intermediates to visualize and debug.\n'
                 '\n'
                 '5. **Improved Pattern Definitions**: The existing patterns '
                 'should be corrected and thoroughly checked by '
                 'cross-referencing with known successful matches.\n'
                 '\n'
                 "Now, let's revise the code to address these issues:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 26,
    'num_false_positives': 100,
    'num_true_negatives': 36512,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.20634920634920634,
    'recall': 1.0,
    'f1': 0.34210526315789475,
    'accuracy': 0.9972705933730007,
    'negative_predictive_value': 1.0}