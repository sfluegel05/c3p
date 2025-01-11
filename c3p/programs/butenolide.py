"""
Classifies: CHEBI:50523 butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide consists of a 5-membered lactone ring (furanone structure), often with a ketone group or 
    related functional group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for butenolide structure
    # More flexible pattern to catch variations in the butenolide structure
    butenolide_pattern_1 = Chem.MolFromSmarts("C1=COC(=O)C1")  # Basic furanone pattern
    butenolide_pattern_2 = Chem.MolFromSmarts("C1=CC(=O)OC1")  # Variation on furanone positioning

    # Check if molecule matches any of the butenolide substructures
    if mol.HasSubstructMatch(butenolide_pattern_1):
        return True, "Contains a basic furanone (butenolide) structure"
    elif mol.HasSubstructMatch(butenolide_pattern_2):
        return True, "Contains a variation of furanone (butenolide) structure"
    
    return False, "Does not contain a recognized butenolide structure"

# Test the function with examples
examples = [
    "CC1(C)OC(=O)C(OCC2CC2)=C1c1ccc(cc1)S(C)(=O)=O",  # firocoxib
    "O=C1OC(=C)C(=C1C(=O)CCCC=CCCCCCC)O",            # Agglomerin B
]

for example in examples:
    result, reason = is_butenolide(example)
    print(f"SMILES: {example} => Result: {result}, Reason: {reason}")


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:50523',
                          'name': 'butenolide',
                          'definition': 'A gamma-lactone that consists of a '
                                        '2-furanone skeleton and its '
                                        'substituted derivatives.',
                          'parents': ['CHEBI:24129', 'CHEBI:37581'],
                          'xrefs': ['Wikipedia:Butenolide'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.0909090909090909 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: O=C1OC(=O)C(=C1/C=C/CCCCCC)C NAME: '
               'Tricladolide D REASON: CORRECT Contains a butenolide '
               '(furanone) structure\n'
               ' * SMILES: O1C(=O)C(CCCCCCCCCCCCCC)=C(C1=O)C NAME: '
               'Chaetomellic acid A anhydride REASON: CORRECT Contains a '
               'butenolide (furanone) structure\n'
               ' * SMILES: O=C1OC(=O)C(=C1/C=C/C=C/C=C(\\C(=O)C)/C)C NAME: '
               'Graphenone REASON: CORRECT Contains a butenolide (furanone) '
               'structure\n'
               ' * SMILES: O1C(=O)C(=C(C1=O)C)C NAME: 2,3-Dimethylmaleic '
               'anhydride REASON: CORRECT Contains a butenolide (furanone) '
               'structure\n'
               ' * SMILES: O1C(=O)C(O)=C(O)C1=O NAME: '
               '3,4-dihydroxy-2,5-furandione REASON: CORRECT Contains a '
               'butenolide (furanone) structure\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'CC1(C)OC(=O)C(OCC2CC2)=C1c1ccc(cc1)S(C)(=O)=O NAME: firocoxib '
               'REASON: MISSED Does not contain a butenolide structure\n'
               ' * SMILES: '
               'O[C@@]12[C@@]([C@](CC1)([C@]3(CC(OC3)=O)[H])[H])([C@H](O)C[C@]4(C2CC[C@]5([C@@]4(CC[C@H](O)C5)C)[H])[H])C '
               'NAME: 20, 22-Dihydrodigoxigenin REASON: MISSED Does not '
               'contain a butenolide structure\n'
               ' * SMILES: '
               'O=C1OC2=C3[C@@]([C@@H](C(C)C)[C@H]([C@H]3O)O)(CC[C@]4(C2=C1[C@@H](O)CC4)C)C '
               'NAME: Plicatilisin A REASON: MISSED Does not contain a '
               'butenolide structure\n'
               ' * SMILES: '
               'CC1C(C(C(C(O1)OC2CC(C3(C4C(CCC3(C2)O)C5(CCC(C5(CC4O)C)C6=CC(=O)OC6)O)CO)O)O)O)O '
               'NAME: '
               '3-[1,5,11,14-tetrahydroxy-10-(hydroxymethyl)-13-methyl-3-[(3,4,5-trihydroxy-6-methyl-2-oxanyl)oxy]-2,3,4,6,7,8,9,11,12,15,16,17-dodecahydro-1H-cyclopenta[a]phenanthren-17-yl]-2H-furan-5-one '
               'REASON: MISSED Does not contain a butenolide structure\n'
               ' * SMILES: O=C1O[C@@H](CC(=O)CC)C(=C1)OC NAME: Annularin H '
               'REASON: MISSED Does not contain a butenolide structure\n'
               ' * SMILES: '
               'O=C1O[C@@H]2[C@@H]3[C@@]([C@@H](C(C)C)[C@H]([C@H]3OC(=O)C)O)(CC[C@]4(C2=C1[C@@H](O)CC4)C)C '
               'NAME: Radianspene I REASON: MISSED Does not contain a '
               'butenolide structure\n'
               ' * SMILES: O=C1O[C@@H]([C@@H](O)C)C=C1[C@H](OC)[C@H](O)C NAME: '
               'Aspilactonol D REASON: MISSED Does not contain a butenolide '
               'structure\n'
               ' * SMILES: '
               'C\\C=C(/C)C(=O)O[C@@H](C[C@@]1(C)[C@H](C)C[C@H](OC(C)=O)[C@]2(COC(C)=O)[C@@H]1CCC[C@]21CO1)C1=CC(=O)OC1 '
               'NAME: ajuganipponin B REASON: MISSED Does not contain a '
               'butenolide structure\n'
               ' * SMILES: O=C1OC(=C)C(=C1C(=O)CCCC=CCCCCCC)O NAME: Agglomerin '
               'B REASON: MISSED Does not contain a butenolide structure\n'
               ' * SMILES: O=C1OC([C@H]2O[C@H](C(O)C)C=C2)C=C1 NAME: '
               'Diplobifuranylone D REASON: MISSED Does not contain a '
               'butenolide structure\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Does not contain a recognized '
                                               'butenolide structure'}],
    'sample_false_negatives': [   {   'smiles': 'O[C@@]12[C@@]([C@](CC1)([C@]3(CC(OC3)=O)[H])[H])([C@H](O)C[C@]4(C2CC[C@]5([C@@]4(CC[C@H](O)C5)C)[H])[H])C',
                                      'name': '20, 22-Dihydrodigoxigenin',
                                      'reason': 'Does not contain a recognized '
                                                'butenolide structure'},
                                  {   'smiles': '[H][C@@]1(OC(=O)C(=O)C1=O)[C@@H](O)CO',
                                      'name': 'L-dehydroascorbic acid',
                                      'reason': 'Does not contain a recognized '
                                                'butenolide structure'},
                                  {   'smiles': '[H][C@@]1(CC[C@]2(O)[C@]3([H])CC[C@]4(O)C[C@H](C[C@@H](O)[C@]4(CO)[C@@]3([H])[C@H](O)C[C@]12C)O[C@@H]1O[C@@H](C)[C@H](O)[C@@H](O)[C@H]1O)C1COC(=O)C1',
                                      'name': 'dihydroouabain',
                                      'reason': 'Does not contain a recognized '
                                                'butenolide structure'},
                                  {   'smiles': 'C(C(C1C(C(C(=O)O1)O)O)=N)O',
                                      'name': '3,4-dihydroxy-5-(2-hydroxyethanimidoyl)oxolan-2-one',
                                      'reason': 'Does not contain a recognized '
                                                'butenolide structure'},
                                  {   'smiles': 'O.[H][C@@]1(OC(=O)C(=O)C1=O)[C@@H](O)CO',
                                      'name': 'L-dehydroascorbic acid hydrate',
                                      'reason': 'Does not contain a recognized '
                                                'butenolide structure'}],
    'attempt': 1,
    'reasoning': "To improve the program for classifying butenolides, let's "
                 'first analyze the issues with the current approach:\n'
                 '\n'
                 '1. **Pattern Specificity**: The SMARTS pattern used for '
                 'detecting butenolides might be too strict or not entirely '
                 'representative of all butenolide structures. Butenolides are '
                 'a subclass of lactones containing a 5-membered ring with a '
                 'certain degree of unsaturation, often represented by a '
                 'furanone-like structure. The structure may vary slightly '
                 'between different butenolides.\n'
                 '\n'
                 '2. **Missed Substructural Diversity**: The lactone ring in '
                 'butenolides can have different substituents around the ring, '
                 'which may not be covered by the current SMARTS pattern. \n'
                 '\n'
                 '3. **Structural Variability**: The current pattern does not '
                 'account for variations such as different positions of double '
                 'bonds or additional groups attached to the core butenolide '
                 'structure.\n'
                 '\n'
                 'To address these issues, we can redefine the pattern to be '
                 'broader while still aiming to capture the essence of a '
                 'butenolide structure. We will incorporate variability in the '
                 'ring and look for the furanone structure with a ketone ring:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 420,
    'num_false_positives': 100,
    'num_true_negatives': 13995,
    'num_false_negatives': 5,
    'num_negatives': None,
    'precision': 0.8076923076923077,
    'recall': 0.9882352941176471,
    'f1': 0.888888888888889,
    'accuracy': 0.9927685950413223,
    'negative_predictive_value': 0.9996428571428572}