"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    A limonoid is a highly oxygenated triterpenoid with a structure commonly
    derived from a 4,4,8-trimethyl-17-furanylsteroid skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for a reasonable range of carbon atoms typical in limonoids
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 40:
        return False, f"Number of carbon atoms {c_count} outside expected range for limonoid"

    # Ensure a high range of oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 4:
        return False, f"Insufficient oxygen content for a highly oxygenated limonoid: {o_count}"

    # Identify oxygen functional groups: ketones, esters, hydroxyl, and potential lactone
    important_oxy_patterns = [
        Chem.MolFromSmarts("C=O"),
        Chem.MolFromSmarts("C(=O)O"),
        Chem.MolFromSmarts("[OX2H]"),  # Alcohol
        Chem.MolFromSmarts("O1CCO1"),  # Broader to catch different heterocycles
        Chem.MolFromSmarts("C1(O)CC1"),  # Another form of closed-ring oxygen groups
    ]
    for oxy_pattern in important_oxy_patterns:
        if mol.HasSubstructMatch(oxy_pattern):
            break
    else:
        return False, "Missing characteristic oxygenated functionalities"

    # Check if there is any form of methyl groups, allowing for variations
    methyl_patterns = [
        Chem.MolFromSmarts("C(C)(C)C"),  # Standard trimethyl
        Chem.MolFromSmarts("C(C)(C)(C)"),  # Various methylated forms
        Chem.MolFromSmarts("C(C)([OX2])")  # Methyl ether or alcohol derivatives
    ]
    if not any(mol.HasSubstructMatch(pat) for pat in methyl_patterns):
        return False, "Lack of characteristic methylation"

    return True, "SMILES string corresponds to a limonoid structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:39434',
                          'name': 'limonoid',
                          'definition': 'Any triterpenoid that is highly '
                                        'oxygenated and has a prototypical '
                                        'structure either containing or '
                                        'derived from a precursor with a '
                                        '4,4,8-trimethyl-17-furanylsteroid '
                                        "skeleton. The term 'limonoid' comes "
                                        'from limonin, the first '
                                        'tetranortriterpenoid obtained from '
                                        'citrus bitter principles.',
                          'parents': ['CHEBI:36615'],
                          'xrefs': ['PMID:16462017', 'Wikipedia:Limonoid'],
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
               'CC(=O)O[C@@H]1C[C@H]2C(C)(C)C(=O)C=C[C@]2(C)[C@H]2CC[C@]3(C)OC(=O)C=C3[C@]12C '
               'NAME: beta-nimolactone REASON: MISSED Lack of furan ring '
               'structure which is common in limonoids\n'
               ' * SMILES: '
               '[H][C@]1(O)OC[C@@]23[C@H](C[C@@H](OC(C)=O)[C@]4(CO[C@@H]([C@@H](O)[C@](C)([C@]12[H])[C@]12O[C@@]1(C)[C@H]1C[C@@H]2O[C@]2([H])OC=C[C@]12O)[C@]34[H])C(=O)OC)OC(=O)C(\\C)=C\\C '
               'NAME: azadirachtin H REASON: MISSED Lack of furan ring '
               'structure which is common in limonoids\n'
               ' * SMILES: '
               'CC(=O)O[C@@H]1C[C@H]2C(C)(C)C(=O)C=C[C@]2(C)[C@H]2CC[C@@]3(C)[C@@H](CC=C3[C@]12C)C1COC(=O)C1 '
               'NAME: 20,21,22,23-tetrahydro-23-oxoazadirone REASON: MISSED '
               'Lack of furan ring structure which is common in limonoids\n'
               ' * SMILES: '
               'COC(=O)C[C@H]1[C@]2(C)C[C@@]3(O)[C@]1(C)[C@H]1CC[C@@]4(C)[C@@H](OC(=O)[C@H](OC(C)=O)C4=C1[C@H](OC(=O)C(C)C)[C@@]3(O)[C@H]2OC(=O)C(\\C)=C\\C)c1ccoc1 '
               'NAME: 15-acetyltrichagmalin C REASON: MISSED Lack of furan '
               'ring structure which is common in limonoids\n'
               ' * SMILES: '
               'CC(=O)O[C@@H]1C[C@H]2C(C)(C)C(=O)C=C[C@]2(C)[C@H]2CC[C@@]3(C)[C@@H](CC=C3[C@]12C)c1ccoc1 '
               'NAME: azadirone REASON: MISSED Insufficient oxygen content for '
               'a highly oxygenated limonoid: 4\n'
               ' * SMILES: '
               '[H][C@@]1(C[C@]([H])(C(O)O1)[C@]1([H])CC[C@]23C[C@]12CC[C@]1([H])[C@@]2(C)CC[C@@H](OC(=O)CC(C)C)C(C)(C)[C@]2([H])C[C@@H](O)[C@@]31C)[C@]1([H])OC1(C)C '
               'NAME: skimmiarepin A REASON: MISSED Lack of furan ring '
               'structure which is common in limonoids\n'
               ' * SMILES: '
               'CO[C@H]1C[C@@]2(O)[C@@H](O1)O[C@H]1C[C@@H]2[C@]2(C)O[C@]12[C@@]1(C)[C@@H]2[C@](O)(OC[C@]22[C@H]3[C@@H](OC[C@@]3([C@@H](C[C@@H]2OC(=O)C(\\C)=C\\C)OC(C)=O)C(=O)OC)[C@H]1O)C(=O)OC '
               'NAME: 23-epivepaol REASON: MISSED Lack of furan ring structure '
               'which is common in limonoids\n'
               ' * SMILES: '
               'COC(=O)C[C@H]1[C@]2(C)C[C@@]3(O)[C@]1(C)[C@H]1CC[C@@]4(C)[C@@H](OC(=O)[C@H](OC(=O)C(C)(C)O)C4=C1[C@H](OC(C)=O)[C@@]3(O)[C@H]2OC(=O)C(\\C)=C\\C)c1ccoc1 '
               'NAME: 30-acetyltrichagmalin F REASON: MISSED Lack of furan '
               'ring structure which is common in limonoids\n'
               ' * SMILES: '
               '[C@]123[C@](C(=O)[C@]([C@@]1(CC[C@@]4([C@]5(C=CC(C([C@@]5(C[C@H]([C@@]24C)OC(=O)C6=CC=CC=C6)[H])(C)C)=O)C)[H])C)(C7=COC=C7)[H])([H])O3 '
               'NAME: 7-deacetyl-7-benzoylepoxyazadiradione REASON: MISSED '
               'Lack of furan ring structure which is common in limonoids\n'
               ' * SMILES: '
               'CO[C@@H]1C[C@@]2(O)[C@@H](O1)O[C@H]1C[C@@H]2[C@]2(C)O[C@]12[C@@]1(C)[C@@H]2[C@](O)(OC[C@]22[C@H]3[C@@H](OC[C@@]3([C@@H](C[C@@H]2OC(=O)C(\\C)=C\\C)OC(C)=O)C(=O)OC)[C@H]1O)C(=O)OC '
               'NAME: vepaol REASON: MISSED Lack of furan ring structure which '
               'is common in limonoids\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'Lack of characteristic '
                                               'methylation'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'Insufficient oxygen content '
                                               'for a highly oxygenated '
                                               'limonoid: 1'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'Number of carbon atoms 9 '
                                               'outside expected range for '
                                               'limonoid'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'Number of carbon atoms 19 '
                                               'outside expected range for '
                                               'limonoid'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'Number of carbon atoms 6 '
                                               'outside expected range for '
                                               'limonoid'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'Number of carbon atoms 42 '
                                               'outside expected range for '
                                               'limonoid'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'Number of carbon atoms 68 '
                                               'outside expected range for '
                                               'limonoid'},
                                 {   'smiles': 'O=C(N(CCO)CCO)CN1C(=NC=C1)[N+]([O-])=O',
                                     'name': '1H-Imidazole-1-acetamide, '
                                             'N,N-bis(2-hydroxyethyl)-2-nitro-',
                                     'reason': 'Number of carbon atoms 9 '
                                               'outside expected range for '
                                               'limonoid'},
                                 {   'smiles': 'CC1=NN(C(C1)(C2=CN=CC=C2)O)C(=O)CC3=CC=C(C=C3)C',
                                     'name': '1-[5-hydroxy-3-methyl-5-(3-pyridinyl)-4H-pyrazol-1-yl]-2-(4-methylphenyl)ethanone',
                                     'reason': 'Number of carbon atoms 18 '
                                               'outside expected range for '
                                               'limonoid'},
                                 {   'smiles': 'COC(=O)C[C@@H]1CC[C@H]([C@H](O1)CO)NC(=O)C2=CC=CC=C2',
                                     'name': '2-[(2S,5R,6S)-5-benzamido-6-(hydroxymethyl)-2-oxanyl]acetic '
                                             'acid methyl ester',
                                     'reason': 'Number of carbon atoms 16 '
                                               'outside expected range for '
                                               'limonoid'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'The initial implementation of the `is_limonoid` function '
                 'failed because it was too broad and rigid. The outcomes '
                 'indicated the following issues:\n'
                 '\n'
                 '1. **Furan Ring Requirement**: The classification required a '
                 'furan ring, expressed as `o1ccc1`, which is common but not '
                 'universal among limonoids. Some examples probably lacked an '
                 'explicit furan ring yet still met the criteria through other '
                 'defining features.\n'
                 '\n'
                 '2. **Trimethyl Group**: The pattern `[C(C)(C)C]` was '
                 'required for classifying a limonoid, though limonoids may '
                 'have diverse modifications to this structural pattern due to '
                 'various substituent rearrangements or oxidation states.\n'
                 '\n'
                 '3. **Oxygen Content and Functionalities**: The high '
                 'oxygenation characteristic needs a nuanced check, ensuring '
                 'not only quantity but variety of oxygen-based '
                 'functionalities were detected which were oversimplified in '
                 'the initial code.\n'
                 '\n'
                 'To improve the classification, the following modifications '
                 'are proposed:\n'
                 '\n'
                 '1. **Remove strict dependency on the explicit furan ring '
                 'presence**, allowing for a broader variety of oxygenated '
                 'heterocyclic rings or highly oxygenated functional '
                 'characteristics.\n'
                 '\n'
                 '2. **Extend the detection of high oxygenation to capture a '
                 'broader range of structural motifs** relevant to limonoids. '
                 'This could include additional SMARTS patterns for esters, '
                 'lactones, and epoxides.\n'
                 '\n'
                 '3. **Use a more flexible requirement for trimethyl groups**, '
                 'acknowledging various potential rearrangements and '
                 'extensions of the carbon framework inherent to limonoids.\n'
                 '\n'
                 '4. **Refine context of carbon and oxygen count ranges**, '
                 'based on limonoid structure paradigms, accommodating '
                 'empirical diversity in these compounds.\n'
                 '\n'
                 "Here's the redesigned function:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 87,
    'num_false_positives': 100,
    'num_true_negatives': 230,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.46524064171123,
    'recall': 1.0,
    'f1': 0.635036496350365,
    'accuracy': 0.7601918465227818,
    'negative_predictive_value': 1.0}