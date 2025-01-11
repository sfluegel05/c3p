"""
Classifies: CHEBI:47916 flavonoid
"""
"""
Classifies: CHEBI:47916 flavonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_flavonoid(smiles: str):
    """
    Determines if a molecule is a flavonoid based on its SMILES string.
    A flavonoid has a 1-benzopyran core with an aryl substituent at position 2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a flavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Core patterns for different flavonoid classes - simplified and more general
    core_patterns = [
        # Basic benzopyran core (most general)
        "O1c2ccccc2CC1",
        
        # Flavone/flavonol core (more general)
        "O1c2ccccc2C(=O)C1",
        
        # Isoflavone core
        "O1c2ccccc2C(=O)C=C1",
        
        # Flavanone core
        "O1CCc2ccccc2C1=O",
        
        # Anthocyanidin core (includes charged form)
        "O1c2ccccc2C=[O+]C1",
        
        # Chalcone-type pattern
        "O=CC(=O)c1ccccc1O",
        
        # More specific but common patterns
        "O1c2c(O)cc(O)cc2OC(c2ccccc2)C1=O",
        "O1c2c(O)cc(O)cc2OC(c2ccccc2)C1",
    ]

    has_core = False
    matched_pattern = None
    for pattern in core_patterns:
        pat = Chem.MolFromSmiles(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            has_core = True
            matched_pattern = pattern
            break

    if not has_core:
        # Try more general SMARTS patterns if SMILES patterns fail
        smarts_patterns = [
            # Very general benzopyran core
            "O1[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]2[#6][#6]1",
            # General chromone core
            "O1[#6]~2~[#6]~[#6]~[#6]~[#6]~[#6]2[#6](=O)[#6]1",
        ]
        for pattern in smarts_patterns:
            pat = Chem.MolFromSmarts(pattern)
            if pat is not None and mol.HasSubstructMatch(pat):
                has_core = True
                matched_pattern = pattern
                break

    if not has_core:
        return False, "No flavonoid core structure found"

    # Check for aromatic/phenyl substituent - more general patterns
    aryl_patterns = [
        "c1ccccc1",  # Basic phenyl
        "c1cccc([#8,#7,#6])c1",  # Substituted phenyl
        "c1cc([#8,#7])cc([#8,#7])c1",  # Di-substituted phenyl
    ]

    has_aryl = False
    for pattern in aryl_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None and mol.HasSubstructMatch(pat):
            has_aryl = True
            break

    # Characteristic flavonoid features
    characteristic_patterns = [
        "[OX2H1]",  # Hydroxyl
        "[OX2][CH3]",  # Methoxy
        "C(=O)",  # Carbonyl
        "O[CH]1O[CH][CH][CH][CH][CH]1",  # Sugar ring
        "O[CH]1O[CH][CH][CH][CH]1",  # Furanose sugar
    ]

    characteristic_count = 0
    for pattern in characteristic_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            matches = len(mol.GetSubstructMatches(pat))
            characteristic_count += matches

    # Ring count check
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 2:
        return False, "Insufficient number of rings"

    # Classification logic
    if has_aryl and characteristic_count >= 1:
        return True, "Contains flavonoid core structure with appropriate substituents"
    elif has_core and characteristic_count >= 2:
        return True, "Contains flavonoid core with characteristic substitution pattern"
    else:
        return False, "Missing required structural features for flavonoid classification"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:47916',
                          'name': 'flavonoid',
                          'definition': "Any member of the 'superclass' "
                                        'flavonoids whose skeleton is based on '
                                        '1-benzopyran with an aryl substituent '
                                        'at position 2. The term was '
                                        'originally restricted to natural '
                                        'products, but is now also used to '
                                        'describe semi-synthetic and fully '
                                        'synthetic compounds.',
                          'parents': ['CHEBI:72544'],
                          'xrefs': ['KEGG:C01579', 'Wikipedia:Flavonoid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.07407407407407407 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O1[C@@H](CC(=O)C=2C1=C(C(O)=C(C2O)CCO)CC=C(C)C)C3=CC(OC)=C(O)C=C3 '
               'NAME: Laxiflorin REASON: CORRECT Contains flavonoid core '
               'structure with appropriate substituents\n'
               ' * SMILES: '
               'O1C2=C(CC(O)C(O)(C)C)C(O)=C(C(O)=C2C(=O)CC1C3=CC=C(O)C=C3)CC=C(C)C '
               'NAME: Lonchocarpol B REASON: CORRECT Contains flavonoid core '
               'structure with appropriate substituents\n'
               ' * SMILES: '
               'CC(C)=CCc1c(O)c(CC=C(C)C)c2O[C@@H](CC(=O)c2c1O)c1cc(O)cc(O)c1 '
               'NAME: '
               '(S)-2-(3,5-dihydroxyphenyl)-5,7-dihydroxy-6,8-bis(3-methylbut-2-enyl)chroman-4-one '
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substituents\n'
               ' * SMILES: '
               'COc1ccc(cc1OC)[C@@H]1CC(=O)c2c(O)c(CC=C(C)C)c3OC(C)(C)C=Cc3c2O1 '
               'NAME: '
               "(2S)-6-(gamma,gamma-dimethylallyl)-5-hydroxy-3',4'-dimethoxy-6'',6''-dimethylpyran[2'',3'':7,8]flavanone "
               'REASON: CORRECT Contains flavonoid core structure with '
               'appropriate substituents\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O1C2=C(C(O)=C(C(O)=C2)C)C(=O)C=C1C3=CC=CC=C3 NAME: '
               'Strobochrysin REASON: MISSED No flavonoid core structure '
               'found\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H](COC(=O)\\C=C\\c2ccc(O)cc2)O[C@@H](Oc2c(oc3cc(O)cc(O)c3c2=O)-c2ccc(O)cc2)[C@H](O)[C@H]1O '
               'NAME: tribuloside REASON: MISSED No flavonoid core structure '
               'found\n'
               ' * SMILES: '
               'O1C(OC=2C=C(C3C4=C(OC(=O)C3)C=5CC(O)C(OC5C=C4O)C6=CC(O)=C(O)C=C6)C=CC2O)C(O)C(O)C(O)C1C(O)=O '
               'NAME: '
               '6-{5-[12-(3,4-dihydroxyphenyl)-8,13-dihydroxy-4-oxo-3,11-dioxatricyclo[8.4.0.0?,?]tetradeca-1(10),2(7),8-trien-6-yl]-2-hydroxyphenoxy}-3,4,5-trihydroxyoxane-2-carboxylic '
               'acid REASON: MISSED No flavonoid core structure found\n'
               ' * SMILES: O=C1C2=C(OC(=C1)C3=C(O)C=C(O)C=C3C)C=C(O)C=C2O '
               'NAME: Penimethavone A REASON: MISSED No flavonoid core '
               'structure found\n'
               ' * SMILES: Oc1cc(O)c2c(c1)oc(cc2=O)-c1cc(O)c(O)cc1O NAME: '
               'isoetin REASON: MISSED No flavonoid core structure found\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](Oc2cc3c(O[C@@H]4O[C@H](COC(=O)\\C=C\\c5ccc(O)cc5)[C@@H](O)[C@H](O)[C@H]4O)cc(O)cc3[o+]c2-c2ccc(O)cc2)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: Pelargonidin 3-O-beta-D-glucoside '
               '5-O-(6-coumaroyl-beta-D-glucoside) REASON: MISSED No flavonoid '
               'core structure found\n'
               ' * SMILES: '
               'O1[C@@H]([C@@H](O)C(OC(=O)/C=C/C2=CC=C(O)C=C2)C(O)[C@@H]1OC3=C(OC=4C(C3=O)=C(O)C=C(O)C4)C5=CC(OC)=C(O)C=C5)COC(=O)/C=C/C6=CC=C(O)C=C6 '
               "NAME: Isorhamnetin 3-(3'',6''-di-p-coumarylglucoside) REASON: "
               'MISSED No flavonoid core structure found\n'
               ' * SMILES: '
               'O([C@H]1[C@@H]([C@H]([C@@H]([C@H](O1)CO[C@]2([C@@H]([C@H]([C@@H]([C@H](O2)CO)O)O)O)[H])O)O)O)C=3C(C4=C(OC3C5=CC=C(O)C=C5)C=C(O[C@@H]6O[C@H]([C@@H]([C@H]([C@H]6O)O)O)C)C=C4O)=O '
               'NAME: kaempferol 3-O-gentiobioside-7-O-rhamnoside REASON: '
               'MISSED No flavonoid core structure found\n'
               ' * SMILES: COc1cc2oc(-c3ccc(O)c(O)c3)c(OC)c(=O)c2c(O)c1O NAME: '
               "3',4',5,6-tetrahydroxy-3,7-dimethoxyflavone REASON: MISSED No "
               'flavonoid core structure found\n'
               ' * SMILES: COc1cc(OC)c2c(c1)oc(-c1ccccc1)c(OC)c2=O NAME: '
               'galangin 3,5,7-trimethyl ether REASON: MISSED No flavonoid '
               'core structure found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O1C(OC=2C=C(C3C4=C(OC(=O)C3)C=5CC(O)C(OC5C=C4O)C6=CC(O)=C(O)C=C6)C=CC2O)C(O)C(O)C(O)C1C(O)=O',
                                      'name': '6-{5-[12-(3,4-dihydroxyphenyl)-8,13-dihydroxy-4-oxo-3,11-dioxatricyclo[8.4.0.0?,?]tetradeca-1(10),2(7),8-trien-6-yl]-2-hydroxyphenoxy}-3,4,5-trihydroxyoxane-2-carboxylic '
                                              'acid',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@H]1O[C@@H](Oc2cc3c(O[C@@H]4O[C@H](COC(=O)\\C=C\\c5ccc(O)cc5)[C@@H](O)[C@H](O)[C@H]4O)cc(O)cc3[o+]c2-c2ccc(O)cc2)[C@H](O)[C@@H](O)[C@@H]1O',
                                      'name': 'Pelargonidin '
                                              '3-O-beta-D-glucoside '
                                              '5-O-(6-coumaroyl-beta-D-glucoside)',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O(C1C(O)C(O)C(OC1OC=2C(=[O+]C=3C(C2)=C(O)C=C(O)C3)C4=CC(O)=C(O)C=C4)COC(=O)C=CC5=CC=C(O)C=C5)C6OC(C(O)C(O)C6O)CO',
                                      'name': '[(6-{[2-(3,4-dihydroxyphenyl)-7-hydroxy-5-oxo-5H-chromen-3-yl]oxy}-3,4-dihydroxy-5-{[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy}oxan-2-yl)methyl][1-hydroxy-3-(4-hydroxyphenyl)prop-2-en-1-ylidene]oxidanium',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'C1(=CC=2O[C@@H]([C@@H](CC2C(=C1)O)OC(C3=CC(=C(C(=C3)O)O)O)=O)C4=CC=C(C=C4)O)O',
                                      'name': 'epiafzelechin 3-gallate',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'CC(C)=CCc1c(O)cc(O)c2C(=O)C[C@H](Oc12)c1ccc(O)cc1O',
                                      'name': 'leachianone G',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O1C(CC(=O)C2=C1C3=C(OC(C=C3)(C)C)C=C2)C4=CC=C(OC)C=C4',
                                      'name': 'Dorspoinsettifolin',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O1[C@@]([C@H](O)CC=2C1=CC(O)=CC2O)(C3=C(O)C=CC(O)=C3)[H]',
                                      'name': '(2R-cis)-2-(2,5-Dihydroxyphenyl)-3,4-dihydro-2H-1-benzopyran-3,5,7-triol',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O1[C@@H](CC(=O)C=2C1=C(C(O)=C(C2O)CCO)CC=C(C)C)C3=CC(OC)=C(O)C=C3',
                                      'name': 'Laxiflorin',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O[C@@H]1Cc2c(O)cc(O)c([C@@H]3[C@@H](O)[C@H](Oc4cc(O)cc(O)c34)c3ccc(O)c(O)c3)c2O[C@@H]1c1ccc(O)c(O)c1',
                                      'name': 'procyanidin B2',
                                      'reason': 'No flavonoid core structure '
                                                'found'},
                                  {   'smiles': 'O=C1C2=C(O)C(=C(O)C=C2O[C@@H](C1)C3=CC=C(OC)C=C3)CC=C(C)C',
                                      'name': '5,7-dihydroxy-2-(4-methoxyphenyl)-6-(3-methylbut-2-enyl)-2,3-dihydrochromen-4-one',
                                      'reason': 'No flavonoid core structure '
                                                'found'}],
    'attempt': 3,
    'reasoning': "Let's analyze what went wrong and how to improve:\n"
                 '\n'
                 '1. Core Issue: The previous code failed to detect many true '
                 'flavonoids, indicating the core structure matching is too '
                 'restrictive. Looking at the false negatives, many are clear '
                 "flavonoids that weren't detected.\n"
                 '\n'
                 '2. Key Problems:\n'
                 '- The SMARTS patterns may be too specific/rigid\n'
                 '- The core pattern matching requires exact matches which may '
                 'fail with complex substitutions\n'
                 '- The characteristic count threshold might be too high\n'
                 '- The bond count threshold of 20 might be excluding some '
                 'valid flavonoids\n'
                 '\n'
                 '3. Improvements needed:\n'
                 '- Simplify and generalize the core flavonoid pattern\n'
                 '- Make the aryl substituent detection more flexible\n'
                 '- Reduce the required characteristic count\n'
                 '- Add more specific patterns for flavonoid subclasses\n'
                 '- Remove the strict bond count requirement\n'
                 '- Add patterns for common modifications like glycosylation\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 94,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.4845360824742268,
    'f1': 0.6527777777777778,
    'accuracy': 0.4845360824742268,
    'negative_predictive_value': 0.0}