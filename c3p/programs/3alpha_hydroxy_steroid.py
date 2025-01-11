"""
Classifies: CHEBI:36835 3alpha-hydroxy steroid
"""
from rdkit import Chem

def is_3alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 3alpha-hydroxy steroid based on its SMILES string.
    A 3alpha-hydroxy steroid is defined as a steroid with a hydroxyl group at the 3-position
    in an alpha orientation.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Flexible steroid backbone pattern, allowing for various stereocenters and saturation.
    steroid_pattern = Chem.MolFromSmarts("C1CCC2C(C(C)CC2)C1")  # Simplified representation of the steroid nucleus
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "No steroid backbone detected"
    
    # Pattern specifically for 3alpha-hydroxy group.
    # The hydroxyl needs to be in the alpha position at the 3-carbon.
    three_alpha_hydroxy_pattern = Chem.MolFromSmarts("[C@@H](O)C(C)C")  # Indicates alpha OH on C3
    if not mol.HasSubstructMatch(three_alpha_hydroxy_pattern):
        return False, "No 3alpha-hydroxy group detected"
    
    return True, "3alpha-hydroxy steroid structure identified"

# Example usage
smiles_example = "[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(CO)CO"
result, reason = is_3alpha_hydroxy_steroid(smiles_example)
print(f"Is 3alpha-hydroxy steroid: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36835',
                          'name': '3alpha-hydroxy steroid',
                          'definition': 'A 3-hydroxy steroid in which the '
                                        '3-hydroxy substituent is in the '
                                        'alpha-position.',
                          'parents': ['CHEBI:35681', 'CHEBI:36834'],
                          'xrefs': [   'MetaCyc:3-alpha-Hydroxysteroids',
                                       'PMID:11514561'],
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
               '[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(CO)CO '
               'NAME: 5beta-cyprinol REASON: MISSED No steroid backbone '
               'detected\n'
               ' * SMILES: '
               '[H][C@@]1(CCC2C3=CC(=O)[C@@]4([H])C[C@H](O)CC[C@]4(C)C3CC[C@]12C)[C@H](C)CCCC(C)COS(O)(=O)=O '
               'NAME: asterasterol B REASON: MISSED No steroid backbone '
               'detected\n'
               ' * SMILES: '
               '[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])CC[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(C)C(O)=O '
               'NAME: 3alpha,7alpha-dihydroxy-5beta-cholestan-26-oic acid '
               'REASON: MISSED No steroid backbone detected\n'
               ' * SMILES: '
               'C1[C@H](C[C@@]2([C@](C1)([C@@]3([C@@]([C@H]([C@H]2O)O)([C@@]4([H])[C@@](CC3)(C)[C@](CC4)([C@@H](CCC(O)=O)C)[H])[H])[H])C)[H])O '
               'NAME: beta-muricholic acid REASON: MISSED No steroid backbone '
               'detected\n'
               ' * SMILES: '
               '[H][C@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCO '
               'NAME: 5alpha-cholane-3alpha,7alpha,12alpha,24-tetrol REASON: '
               'MISSED No steroid backbone detected\n'
               ' * SMILES: '
               '[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CC[C@H](O)C(C)C '
               'NAME: (24S)-5beta-cholestane-3alpha,7alpha,12alpha,24-tetrol '
               'REASON: MISSED No steroid backbone detected\n'
               ' * SMILES: '
               'S(O)(=O)(=O)CCNC(=O)CC[C@H]([C@@]1([C@@]2([C@]([C@]3([C@@]([C@@]4([C@]([C@H](O)[C@H]3O)(C[C@H](O)CC4)[H])C)(CC2)[H])[H])(CC1)[H])C)[H])C '
               'NAME: tauro-alpha-muricholic acid REASON: MISSED No steroid '
               'backbone detected\n'
               ' * SMILES: '
               '[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(CO)COS(O)(=O)=O '
               'NAME: 5beta-cyprinol sulfate REASON: MISSED No steroid '
               'backbone detected\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])COC(=O)[C@@]4([H])C[C@H](O)[C@H](O)C[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)[C@@H](O)[C@H](O)[C@@H](C)C(C)C '
               'NAME: brassinolide REASON: MISSED No steroid backbone '
               'detected\n'
               ' * SMILES: '
               'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3([C@H]([C@H]([C@@]2(C[C@@H](C1)O)[H])O)O)[H])(CC[C@@]4([C@@H](CCC(NCCS(O)(=O)=O)=O)C)[H])[H])C)[H])C '
               'NAME: tauro-beta-muricholic acid REASON: MISSED No steroid '
               'backbone detected\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No steroid backbone detected'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No steroid backbone detected'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]1(CCC2C3=CC(=O)[C@@]4([H])C[C@H](O)CC[C@]4(C)C3CC[C@]12C)[C@H](C)CCCC(C)COS(O)(=O)=O',
                                      'name': 'asterasterol B',
                                      'reason': 'No 3alpha-hydroxy group '
                                                'detected'},
                                  {   'smiles': '[H][C@@]12CC[C@]3([H])[C@]([H])(CC[C@]4(C)[C@@H](OC(C)=O)[C@H](C[C@@]34[H])[N+]3(CCCC3)CC=C)[C@@]1(C)C[C@@H]([C@@H](O)C2)N1CCOCC1',
                                      'name': 'rocuronium',
                                      'reason': 'No steroid backbone detected'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC[C@]4([H])C[C@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)C(C)=O',
                                      'name': '3alpha-hydroxy-5beta-pregnan-20-one',
                                      'reason': 'No 3alpha-hydroxy group '
                                                'detected'},
                                  {   'smiles': 'O1[C@@H]([C@H]([C@@H]([C@H]([C@@H]1O[C@@H]2[C@@]3([C@]([C@]4([C@](CC3)([C@@]5([C@](C[C@@H](CC5)O)(CC4)[H])C)[H])[H])(CC2)[H])C)O)O)O)C(O)=O',
                                      'name': '5alpha-androstane-3alpha,17beta-diol '
                                              '17-glucosiduronic acid',
                                      'reason': 'No steroid backbone detected'},
                                  {   'smiles': '[H][C@@]1(CCC2C3=CC(=O)[C@@]4([H])C[C@H](O)CC[C@]4(C)C3CC[C@]12C)[C@H](C)\\C=C\\CC(C)COS(O)(=O)=O',
                                      'name': 'asterasterol C',
                                      'reason': 'No 3alpha-hydroxy group '
                                                'detected'},
                                  {   'smiles': 'C1[C@@]2([C@@]([C@@]3([C@](C[C@H](O)CC3)(C1)[H])C)(CC[C@@]4([C@H](CC[C@@]24[H])OS(O)(=O)=O)C)[H])[H]',
                                      'name': '(3alpha,5alpha,17beta)-3-hydroxyandrostan-17-yl '
                                              'sulfate',
                                      'reason': 'No steroid backbone detected'},
                                  {   'smiles': 'C[C@H](CCC(O)=O)[C@H]1CC[C@H]2[C@H]3[C@H](CC[C@]12C)[C@@]1(C)CC[C@@H](O)C[C@H]1CC3=O',
                                      'name': '7-oxolithocholic acid',
                                      'reason': 'No 3alpha-hydroxy group '
                                                'detected'},
                                  {   'smiles': 'C[C@H](CCC[C@H](C)C(O)=O)[C@H]1CC[C@H]2C3=CC[C@H]4C[C@H](O)CC[C@]4(C)[C@H]3CC[C@]12C',
                                      'name': '(25S)-3alpha-hydroxy-5alpha-cholest-7-en-26-oic '
                                              'acid',
                                      'reason': 'No 3alpha-hydroxy group '
                                                'detected'},
                                  {   'smiles': 'O[C@H]1C[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@](CC4)([C@@H](CCCC(O)=O)C)[H])(CC3)C)[H])(CC2)[H])[H])(CC1)C)[H]',
                                      'name': 'Homolithocholic acid',
                                      'reason': 'No 3alpha-hydroxy group '
                                                'detected'},
                                  {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C(C[C@@]2(C[C@@H](C1)O)[H])=O)[H])(CC[C@@]4([C@@H](CCC(NCCS(O)(=O)=O)=O)C)[H])[H])C)[H])C',
                                      'name': '7-oxotaurolithocholic acid',
                                      'reason': 'No 3alpha-hydroxy group '
                                                'detected'}],
    'attempt': 4,
    'reasoning': 'The previous attempt to classify 3alpha-hydroxy steroids '
                 'failed because it did not adequately capture the structural '
                 'features that define this class of compounds. Specifically, '
                 'the steroid backbone pattern was too rigid and not broadly '
                 'defined enough to account for the diversity in steroid '
                 'structures. Moreover, the stereo-specific requirement for '
                 'the hydroxyl group at the 3-position was potentially too '
                 'simplistic, neglecting variant ring junctions and '
                 'substitutions naturally occurring in steroids.\n'
                 '\n'
                 "Here's what we can do to improve the classification:\n"
                 '\n'
                 '1. **Broadening the Steroid Backbone Definition**: Instead '
                 'of using a single, rigid SMARTS pattern to identify the '
                 'steroidal nucleus, we should use a more flexible pattern '
                 'that accounts for typical variations in steroid structures, '
                 'such as substitutions and ring junctions.\n'
                 '\n'
                 '2. **Improving Hydroxy Group Recognition**: We need to '
                 'ensure that we identify the 3alpha-hydroxy group by '
                 'explicitly anchoring the alcohol at the 3-position while '
                 'considering common ring-indexing in steroids.\n'
                 '\n'
                 "3. **Adding Stereocenter Verification**: Ensure that we're "
                 'correctly determining the alpha positioning of the hydroxy '
                 'group by checking the specific 3-position stereochemical '
                 'orientation.\n'
                 '\n'
                 "Here's a refined version of the program using these "
                 'improvements:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 72,
    'num_false_positives': 100,
    'num_true_negatives': 4776,
    'num_false_negatives': 19,
    'num_negatives': None,
    'precision': 0.4186046511627907,
    'recall': 0.7912087912087912,
    'f1': 0.5475285171102662,
    'accuracy': 0.9760418763841353,
    'negative_predictive_value': 0.9960375391032326}