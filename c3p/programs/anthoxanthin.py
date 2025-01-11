"""
Classifies: CHEBI:192499 anthoxanthin
"""
"""
Classifies: anthoxanthin
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_anthoxanthin(smiles: str):
    """
    Determines if a molecule is an anthoxanthin based on its SMILES string.
    Anthoxanthins are a type of flavonoid pigments in plants, including flavones and flavonols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an anthoxanthin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Aromaticity perception (ensure correct handling of aromatic rings)
    Chem.SanitizeMol(mol)

    # Define SMARTS patterns for flavone and flavonol cores
    # Flavone core: 2-phenylchromen-4-one
    flavone_pattern = Chem.MolFromSmarts('c1cc(-c2oc3ccccc3c(=O)c2)ccc1')

    # Flavonol core: 3-hydroxyflavone
    flavonol_pattern = Chem.MolFromSmarts('c1cc(-c2oc3ccccc3c(=O)c2O)ccc1')

    # Check for flavone core
    if mol.HasSubstructMatch(flavone_pattern):
        return True, "Contains flavone core structure"

    # Check for flavonol core
    if mol.HasSubstructMatch(flavonol_pattern):
        return True, "Contains flavonol core structure"

    # Define patterns for glycosylated flavones/flavonols (O-glycosides)
    # Simplified pattern: flavone or flavonol core with any sugar attached via oxygen
    o_glycoside_pattern = Chem.MolFromSmarts('[$([cO][CX4H]),$([cO][CX4H][CX4H]),$([cO][CX4H][CX4H][CX4H])]')
    flavone_oglycoside_pattern = Chem.CombineMols(flavone_pattern, o_glycoside_pattern)
    flavonol_oglycoside_pattern = Chem.CombineMols(flavonol_pattern, o_glycoside_pattern)

    if mol.HasSubstructMatch(flavone_oglycoside_pattern):
        return True, "Contains flavone O-glycoside structure"

    if mol.HasSubstructMatch(flavonol_oglycoside_pattern):
        return True, "Contains flavonol O-glycoside structure"

    # Check for C-glycosides (sugar attached directly to carbon)
    # C-glycoside pattern: flavone or flavonol core with sugar attached to C6 or C8
    c_glycoside_pattern = Chem.MolFromSmarts('c1cc(-c2oc3ccccc3c(=O)c2[c,C])ccc1')
    if mol.HasSubstructMatch(c_glycoside_pattern):
        return True, "Contains flavone or flavonol C-glycoside structure"

    # Check for methylated derivatives (methoxy groups)
    methoxy_pattern = Chem.MolFromSmarts('c-oc')
    if mol.HasSubstructMatch(flavone_pattern) or mol.HasSubstructMatch(flavonol_pattern):
        if mol.HasSubstructMatch(methoxy_pattern):
            return True, "Contains methylated flavone or flavonol core"

    # Check for sulfonated derivatives (sulfo groups)
    sulfo_pattern = Chem.MolFromSmarts('S(=O)(=O)[O-]')
    if mol.HasSubstructMatch(sulfo_pattern):
        if mol.HasSubstructMatch(flavone_pattern) or mol.HasSubstructMatch(flavonol_pattern):
            return True, "Contains sulfonated flavone or flavonol"

    # No anthoxanthin core structure found
    return False, "Does not contain anthoxanthin core structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:192499',
                          'name': 'anthoxanthin',
                          'definition': 'are a type of flavonoid pigments in '
                                        'plants. Anthoxanthins are '
                                        'water-soluble pigments which range in '
                                        'color from white or colorless to a '
                                        'creamy to yellow, often on petals of '
                                        'flowers.',
                          'parents': ['CHEBI:47916'],
                          'xrefs': [   'Wikipedia:https://en.wikipedia.org/wiki/Anthoxanthin'],
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
               'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OC=2C(=O)C=3C(OC2C4=CC=C(O)C(=C4)O)=CC(O[C@H]5O[C@H]([C@H](O)[C@@H](O)[C@H]5O)C)=CC3O)CO '
               'NAME: quercetin '
               '3-O-beta-D-glucopyranosyl-7-O-alpha-L-rhamnopyranoside REASON: '
               'MISSED Does not contain anthoxanthin core structure\n'
               ' * SMILES: '
               'OC[C@@H]1O[C@H](Oc2c(oc3cc(O)cc(O)c3c2=O)-c2ccc(O)cc2)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: kaempferol 3-O-beta-L-glucopyranoside REASON: MISSED '
               'Does not contain anthoxanthin core structure\n'
               ' * SMILES: '
               'O[C@H]1[C@H](Oc2cc(O)cc3oc(cc(=O)c23)-c2ccc(O)c(O)c2)O[C@@H]([C@@H](O)[C@@H]1O)C(O)=O '
               'NAME: luteolin 5-O-glucuronide REASON: MISSED Does not contain '
               'anthoxanthin core structure\n'
               ' * SMILES: COc1cc2oc(cc(=O)c2c(OC)c1OC)-c1ccc(O)cc1 NAME: '
               '2-(4-hydroxyphenyl)-5,6,7-trimethoxy-4H-1-benzopyran-4-one '
               'REASON: MISSED Does not contain anthoxanthin core structure\n'
               ' * SMILES: COc1c(OC)c(OC)c2c(oc(cc2=O)-c2ccc(OC(C)=O)cc2)c1OC '
               'NAME: '
               '2-[4-(acetyloxy)phenyl]-5,6,7,8-tetramethoxy-4H-1-benzopyran-4-one '
               'REASON: MISSED Does not contain anthoxanthin core structure\n'
               ' * SMILES: '
               'OC[C@H]1O[C@@H](Oc2cc(O)cc3oc(-c4ccc(O)cc4)c(O)c(=O)c23)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: kaempferol 5-O-beta-D-glucopyranoside REASON: MISSED '
               'Does not contain anthoxanthin core structure\n'
               ' * SMILES: '
               'CC(=O)OC1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C=C3O2)OC(=O)C)OC(=O)C '
               'NAME: acetic acid '
               '[4-(5,7-diacetyloxy-4-oxo-1-benzopyran-2-yl)phenyl] ester '
               'REASON: MISSED Does not contain anthoxanthin core structure\n'
               ' * SMILES: '
               'O1C=2C(=C(O)C=C(OC(=O)C3=CC=CC=C3)C2)C(=O)C=C1C4=CC=CC=C4 '
               'NAME: 5,7-Dihydroxyflavone 7-benzoate REASON: MISSED Does not '
               'contain anthoxanthin core structure\n'
               ' * SMILES: '
               'C=12C(=C(C=C(C1C(C=C(O2)C3=CC=C(C(=C3)O)O)=O)O)O)[C@H]4[C@@H]([C@H]([C@@H]([C@H](O4)CO)O)O)OC5O[C@H]([C@@H]([C@H]([C@H]5O)O)O)C '
               "NAME: orientin 2''-O-rhamnoside REASON: MISSED Does not "
               'contain anthoxanthin core structure\n'
               ' * SMILES: '
               'OC[C@H]1O[C@H]([C@H](O)[C@@H](O)[C@@H]1O)c1c(O)cc2oc(cc(=O)c2c1O)-c1ccc(O)cc1 '
               'NAME: isovitexin REASON: MISSED Does not contain anthoxanthin '
               'core structure\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=C(C=C3)C(=O)N(C)C)O[C@@H]1CN(C)C(=O)NC4=CC=CC=C4F)[C@H](C)CO',
                                     'name': '4-[(4R,5S)-5-[[[(2-fluoroanilino)-oxomethyl]-methylamino]methyl]-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-8-yl]-N,N-dimethylbenzamide',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': 'C[C@H]1[C@H](CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@H]3[C@H](CC[C@@H](O3)CCN(C1=O)C)OC)C)OC',
                                     'name': 'LSM-10564',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': 'CCC(C(=O)NC1=NN=C(S1)SCC)(C(F)(F)F)C(F)(F)F',
                                     'name': 'N-[5-(ethylthio)-1,3,4-thiadiazol-2-yl]-2,2-bis(trifluoromethyl)butanamide',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': 'O=C(N[C@@H](CCC(=O)N)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC1=CC=C(O)C=C1)CC(=O)N',
                                     'name': 'Tyr-Asn-Gln',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': 'OC[C@H]1OC(O)[C@@H](O)[C@@H]1O',
                                     'name': 'D-arabinofuranose',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': 'O[C@@H]([C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)CCC3=C)[H])C)[H])C)C',
                                     'name': '(5Z,7E)-(3S,22R)- '
                                             '24-nor-9,10-seco-5,7,10(19)-cholatriene-3,22-diol',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': 'O1C(C(CC1=O)C)CCCCCC',
                                     'name': '5-Hexyldihydro-4-methyl-2(3H)-furanone',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': 'O=C1O[C@]23C(=O)N[C@H]([C@@H]2[C@H](C)C(=C[C@@H]3C=CC[C@@H]([C@H]([C@H](C=C1)C)O)C)C)CC4=CC=CC=C4',
                                     'name': '7-deoxy-cytochalasin Z7',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': '[C@]123[C@H](O1)[C@](C(C=C2[C@]4([C@](CC3)(C([C@@H](C(C4)=O)O)(C)C)[H])C)=O)(C=C)C',
                                     'name': 'yucalexin P8',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'},
                                 {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(O)=O)[C@@]4(C)[C@@H](O)C[C@]3([H])[C@@]1(C)CC[C@@H](O)C2',
                                     'name': 'deoxycholic acid',
                                     'reason': 'Does not contain anthoxanthin '
                                               'core structure'}],
    'sample_false_negatives': [   {   'smiles': 'Oc1ccc(cc1)[C@H]1Oc2cc(O)cc(O)c2C(=O)[C@@H]1[C@H]1[C@@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)cc1',
                                      'name': 'isochamaejasmin',
                                      'reason': 'Does not contain anthoxanthin '
                                                'core structure'},
                                  {   'smiles': 'C1(OC2=CC(O)=C3C(OC=C(C3=O)C4=CC=C(O)C=C4)=C2)OC(COC(C)=O)C(O)C(O)C1O',
                                      'name': "6''-O-acetylgenistin",
                                      'reason': 'Does not contain anthoxanthin '
                                                'core structure'},
                                  {   'smiles': 'O[C@@H]1[C@@H](O)[C@@H](O[C@@H]([C@H]1O)C(O)=O)Oc1cc(O)c2C(=O)C[C@@H](Oc2c1)c1ccc(O)c(O)c1',
                                      'name': '(2R)-eriodictoyl-7-O-beta-D-glucopyranosiduronic '
                                              'acid',
                                      'reason': 'Does not contain anthoxanthin '
                                                'core structure'},
                                  {   'smiles': 'OC[C@@]1(O)CO[C@@H](OC[C@H]2O[C@H]([C@H](O)[C@@H](O)[C@@H]2O)c2c(O)ccc3c2occ(-c2ccc(O)cc2)c3=O)[C@@H]1O',
                                      'name': 'mirificin',
                                      'reason': 'Does not contain anthoxanthin '
                                                'core structure'},
                                  {   'smiles': 'COC1=CC2OC(=CC(=O)C2C(O)=C1[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)c1ccc(O)c(O)c1',
                                      'name': 'swertiajaponin',
                                      'reason': 'Does not contain anthoxanthin '
                                                'core structure'},
                                  {   'smiles': 'C[C@@H]1O[C@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2Oc2ccc3C(=O)CC(Oc3c2)c2ccccc2)[C@H](O)[C@H](O)[C@H]1O',
                                      'name': '7-hydroxyflavanone '
                                              '7-O-[beta-L-rhamnosyl-(1->2)-beta-D-glucoside]',
                                      'reason': 'Does not contain anthoxanthin '
                                                'core structure'},
                                  {   'smiles': 'OC[C@H]1O[C@@H](Oc2cc(O)c3C(=O)C[C@H](Oc3c2[C@@H]2OC[C@H](O)[C@H](O)[C@H]2O)c2ccc(O)cc2)[C@H](O)[C@@H](O)[C@@H]1O',
                                      'name': '8-C-alpha-L-arabinopyranosyl-7-O-beta-D-glucopyranosylnaringenin',
                                      'reason': 'Does not contain anthoxanthin '
                                                'core structure'},
                                  {   'smiles': '[C@H]1(C=2C3=C(C(C(C4=CC=C(O)C=C4)=CO3)=O)C=CC2[O-])[C@@H]([C@H]([C@H](O)[C@@H](CO)O1)O)O',
                                      'name': 'puerarin(1-)',
                                      'reason': 'Does not contain anthoxanthin '
                                                'core structure'}],
    'attempt': 1,
    'reasoning': 'Upon analyzing the previous code and the outcomes, it '
                 'appears that the SMARTS patterns used to identify '
                 'anthoxanthins were too narrow and did not account for the '
                 'structural diversity of anthoxanthins. Anthoxanthins '
                 'encompass a broad class of flavonoid pigments, including '
                 'flavones, flavonols, and their derivatives such as '
                 'glycosides, methoxylated, and sulfonated forms.\n'
                 '\n'
                 'The previous code only checks for specific flavone and '
                 'flavonol core structures, which may not match the actual '
                 'structures of the provided molecules. Additionally, the '
                 'patterns may not have properly accounted for aromaticity and '
                 'substituents in RDKit.\n'
                 '\n'
                 'To improve the classification, we need to:\n'
                 '\n'
                 '1. **Define more comprehensive SMARTS patterns** that '
                 'accurately represent the flavone and flavonol cores, '
                 'including variations in substituents like hydroxyl and '
                 'methoxy groups.\n'
                 '2. **Account for glycosylated forms** by allowing sugars '
                 'attached at various positions on the core structure.\n'
                 '3. **Ensure correct representation of aromaticity** in the '
                 "SMARTS patterns to match RDKit's handling of aromatic atoms "
                 'and bonds.\n'
                 '4. **Include patterns for other anthoxanthin subclasses**, '
                 'such as flavanones and chalcones, if necessary.\n'
                 '\n'
                 'By broadening the SMARTS patterns and making them more '
                 'flexible, we can improve the detection of anthoxanthins in '
                 'the provided examples.\n'
                 '\n'
                 "Now, here's the updated code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 509,
    'num_false_positives': 100,
    'num_true_negatives': 12102,
    'num_false_negatives': 8,
    'num_negatives': None,
    'precision': 0.8357963875205254,
    'recall': 0.9845261121856866,
    'f1': 0.9040852575488456,
    'accuracy': 0.9915087664124538,
    'negative_predictive_value': 0.9993393889347647}