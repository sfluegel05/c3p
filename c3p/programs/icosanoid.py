"""
Classifies: CHEBI:23899 icosanoid
"""
"""
Classifies: CHEBI:24913 icosanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_icosanoid(smiles: str):
    """
    Determines if a molecule is an icosanoid based on its SMILES string.
    Icosanoids are signaling molecules derived from C20 essential fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an icosanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count oxygens (should have multiple from oxidation)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Too few oxygen atoms for an icosanoid"

    # Look for key structural features
    features_found = []
    
    # Carbon chain patterns characteristic of icosanoids
    chain_patterns = [
        'C~C~C~C~C~C~C~C~C~C~C~C', # Basic long chain
        'C~C=C~C=C~C~C', # Conjugated double bonds
        'C1CCCC1C~C~C~C', # Prostaglandin-like
    ]
    
    for pattern in chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            features_found.append("characteristic carbon chain")
            break

    if not features_found:
        return False, "No characteristic icosanoid carbon chain found"

    # Oxidation and functional group patterns
    patterns = {
        'C=C': "double bonds",
        '[OH]': "hydroxyl groups",
        'C1CCCC1': "cyclopentane ring",
        'C1OC1': "epoxide group",
        'C(=O)C': "ketone group",
        'C(=O)[OH]': "carboxylic acid",
        'C(=O)O[CH2,CH3]': "ester group",
        '[OO]': "peroxide group",
        'C(=O)N': "amide group",
    }
    
    for smart, feature in patterns.items():
        if mol.HasSubstructMatch(Chem.MolFromSmarts(smart)):
            features_found.append(feature)

    # Need multiple oxidation features to be an icosanoid
    oxidation_features = ["hydroxyl groups", "epoxide group", "ketone group", 
                         "carboxylic acid", "peroxide group"]
    oxidation_count = sum(1 for f in oxidation_features if f in features_found)
    
    if oxidation_count < 1:
        return False, "Lacks typical icosanoid oxidation features"

    # Look for specific icosanoid substructures
    icosanoid_patterns = [
        'C(=O)CCC=CCC1C(O)CC(O)C1', # Prostaglandin core
        'C=CC=CC=CC(O)', # Characteristic hydroxylated polyene
        'CC=CC=CC=CC(OO)', # Hydroperoxy pattern
        'C1OC1CC=CC=C', # Epoxy-fatty acid pattern
    ]
    
    for pattern in icosanoid_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            features_found.append("characteristic icosanoid substructure")
            break

    # If we have enough characteristic features, classify as icosanoid
    if len(set(features_found)) >= 3:  # Need at least 3 different features
        features_str = ", ".join(set(features_found))
        return True, f"Contains characteristic icosanoid features: {features_str}"
    
    return False, "Insufficient characteristic icosanoid features"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23899',
                          'name': 'icosanoid',
                          'definition': 'Any member of the group of signalling '
                                        'molecules arising from oxidation of '
                                        'the three C20 essential fatty acids '
                                        '(EFAs) icosapentaenoic acid (EPA), '
                                        'arachidonic acid (AA) and '
                                        'dihomo-gamma-linolenic acid (DGLA).',
                          'parents': ['CHEBI:61697'],
                          'xrefs': [   'LIPID_MAPS_class:LMFA03',
                                       'PMID:12697726',
                                       'PMID:16107763',
                                       'PMID:19160658',
                                       'PMID:2657063',
                                       'Wikipedia:Eicosanoid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7395833333333334 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O[C@H]1[C@@H]([C@@H]([C@@H](O)C1)/C=C/[C@@H](O)C/C=C\\C/C=C\\CC)C/C=C\\CCC(O)=O '
               'NAME: ent-14-epi-14-F4t-NeuroP REASON: CORRECT Contains '
               'characteristic icosanoid features: double bonds, hydroxyl '
               'groups, cyclopentane ring, ketone group\n'
               ' * SMILES: C(\\CC(/C=C/C=C\\C/C=C\\CCCCCO)O)=C\\CCCC(O)=O '
               'NAME: 8,20-DiHETE REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: '
               'CCCCC[C@H](O)\\C=C\\[C@H]1OC(O)C[C@H](O)[C@@H]1CCC(O)=O NAME: '
               '2,3,4,5-tetranor-thromboxane B1 REASON: CORRECT Contains '
               'characteristic icosanoid features: double bonds, hydroxyl '
               'groups, ketone group\n'
               ' * SMILES: '
               'O[C@H]1[C@H]([C@@H](CCCCCCCC(O)=O)C(=O)C1)/C=C/[C@@H](O)CC '
               'NAME: 16-E1c-PhytoP REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, '
               'cyclopentane ring, ketone group\n'
               ' * SMILES: O1C(C(O)CC1C(O)/C=C/C(O)CC)CCCCCCCC(O)=O NAME: '
               'delta14-9-PhytoF REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: '
               'O=C1[C@@H]([C@@H](C/C=C\\CCCCCC(O)=O)C=C1)/C=C/[C@@H](O)CCCCC '
               'NAME: 1a,1b-Dihomo-PGJ2 REASON: CORRECT Contains '
               'characteristic icosanoid features: double bonds, hydroxyl '
               'groups, ketone group\n'
               ' * SMILES: C(CCCO)C/C=C\\CC1C(C/C=C\\C/C=C\\CCCC(O)=O)O1 NAME: '
               '11,12-epoxy-20-hydroxy-(5Z,8Z,14Z)-icosatrienoic acid REASON: '
               'CORRECT Contains characteristic icosanoid features: double '
               'bonds, hydroxyl groups, epoxide group, ketone group\n'
               ' * SMILES: '
               'CCCCC[C@H](O)\\C=C\\[C@@H](C=O)[C@@H](C\\C=C/CCCC(O)=O)C(C)=O '
               'NAME: levuglandin E2 REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: '
               'CCCCC[C@H](O)\\C=C\\C=C/C=C/C=C/[C@H](O)[C@@H](O)CCCC(O)=O '
               'NAME: 6-epi-lipoxin A4 REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: CCCCCCCCC#CCC#CCC#CCCCC(O)=O NAME: '
               '5,8,11-icosatriynoic acid REASON: CORRECT Contains '
               'characteristic icosanoid features: hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: '
               'O[C@H]1[C@@H]([C@@H](C(=O)C1)/C=C/[C@@H](O)CCCCCCCC(O)=O)CC '
               'NAME: 9-E1t-PhytoP REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, '
               'cyclopentane ring, ketone group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/[C@H](O)[C@H]1O[C@H]1C\\C=C/CCCC(O)=O NAME: '
               '(8S,9S)-epoxy-(10R)-hydroxyicosa-(5Z,11Z,14Z)-trienoic acid '
               'REASON: CORRECT Contains characteristic icosanoid features: '
               'double bonds, hydroxyl groups, epoxide group, ketone group\n'
               ' * SMILES: O=C1[C@H]([C@H](CC1)CCCCCCCC(O)=O)CC=CCC NAME: '
               'Oxo-(pentenyl)cyclopentaneoctanoic acid REASON: CORRECT '
               'Contains characteristic icosanoid features: double bonds, '
               'hydroxyl groups, cyclopentane ring, ketone group\n'
               ' * SMILES: '
               'C(\\[C@H](CCCC(O)=O)O)=C\\C=C\\C=C\\[C@@H](C\\C=C/C=C/[C@@H](CC)O)O '
               'NAME: resolvin E1 REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: '
               'C1(=C(C(CC1)=O)C/C=C\\CCCC(O)=O)/C=C/[C@H](CCCC(C)O)O NAME: '
               '19-hydroxy-PGB2 REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: '
               '[C@H](CCCC(O)=O)(/C=C/C=C\\C[C@H](/C=C/C=C\\CCCCC)O)O NAME: '
               '5(S),11(R)-DiHETE REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: '
               'O[C@@H](C(CCCC)(C)C)\\C=C\\[C@@H]1[C@H](C(=O)CC1)C/C=C\\CCCC(O)=O '
               'NAME: 11-Deoxy-16,16-dimethyl-PGE2 REASON: CORRECT Contains '
               'characteristic icosanoid features: double bonds, hydroxyl '
               'groups, cyclopentane ring, ketone group\n'
               ' * SMILES: '
               'O=C1[C@@H]([C@H](C/C=C\\C/C=C\\C/C=C\\CC)C=C1)/C=C/[C@H](O)CCC(O)=O '
               'NAME: 4-epi-4-A4-NeuroP REASON: CORRECT Contains '
               'characteristic icosanoid features: double bonds, hydroxyl '
               'groups, ketone group\n'
               ' * SMILES: '
               'O[C@H]1[C@@H]([C@@H](CCCCCCC(O)=O)C(=O)C1)/C=C/[C@@H](O)CCC[C@H](O)C '
               'NAME: 19R-Hydroxy-PGE1 REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, '
               'cyclopentane ring, ketone group\n'
               ' * SMILES: CCCCC\\C=C/CC1OC1C\\C=C/C\\C=C/CCCC(O)=O NAME: '
               '11,12-EET REASON: CORRECT Contains characteristic icosanoid '
               'features: double bonds, hydroxyl groups, epoxide group, ketone '
               'group\n'
               ' * SMILES: '
               '[C@H](CCCC(O)=O)(/C=C/C=C\\C/C=C\\C\\C=C/C=C/[C@@H](CC)O)OO '
               'NAME: '
               '5(S)-hydroperoxy-18(R)-hydroxy-(6E,8Z,11Z,14Z,16E)-icosapentaenoic '
               'acid REASON: CORRECT Contains characteristic icosanoid '
               'features: double bonds, hydroxyl groups, ketone group\n'
               ' * SMILES: CCCCC\\C=C/C\\C=C/C=C/[C@@H](C\\C=C/CCCC(O)=O)OO '
               'NAME: 8(R)-HPETE REASON: CORRECT Contains characteristic '
               'icosanoid features: double bonds, hydroxyl groups, ketone '
               'group\n'
               ' * SMILES: '
               'C([C@H](/C=C/C=C/C=C\\[C@H](CCCC(O)=O)O)O)/C=C\\CCCCC(O)O '
               'NAME: 20,20-dihydroxyleukotriene B4 REASON: CORRECT Contains '
               'characteristic icosanoid features: double bonds, hydroxyl '
               'groups, ketone group\n'
               ' * SMILES: O[C@H](C/C=C\\CCCCC)/C=C/C=C/C=C/C(=O)CCCC(O)=O '
               'NAME: 5-Oxo-6-trans-leukotriene B4 REASON: CORRECT Contains '
               'characteristic icosanoid features: double bonds, hydroxyl '
               'groups, ketone group\n'
               ' * SMILES: '
               'O=C1C([C@@]2([C@@](C(C(=O)CC2)CCCC)(C1)[H])[H])CCCCCCC(O)=O '
               'NAME: Bicyclo Prostaglandin E1 REASON: CORRECT Contains '
               'characteristic icosanoid features: hydroxyl groups, '
               'cyclopentane ring, ketone group\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O=C1[C@H]([C@@]2([C@H]([C@](C(=O)O)([C@@H](O)CC2)C)CC1)C)C[C@H](OC(=O)C)[C@@]3(C(=O)CC[C@H]3[C@@H](CCC(=C)C(C)C)C)C '
               'NAME: Gloeophyllin J REASON: MISSED Carbon count (32) outside '
               'typical range for icosanoids\n'
               ' * SMILES: '
               'O(C(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)C[C@H](O)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\CCCCC '
               'NAME: DG(18:3n6/0:0/20:4n3) REASON: MISSED Carbon count (41) '
               'outside typical range for icosanoids\n'
               ' * SMILES: '
               'P(OC[C@H](OC(=O)CCC\\C=C\\C[C@@H]1[C@H]([C@H](O)C[C@@H]1O)/C=C/[C@@H](O)CCCCC)COC(=O)CCCCCCCCCCCCCCCCCCCCCCC)(OC[C@H](N)C(O)=O)(O)=O '
               'NAME: PS(24:0/PGF2alpha) REASON: MISSED Carbon count (50) '
               'outside typical range for icosanoids\n'
               ' * SMILES: '
               'O([C@]1(C/C=C\\CCCCC)/C(/C(=O)C=C1)=C/C=C/[C@@H]2OC(=O)CC2)C(=O)C '
               'NAME: clavulolactone III REASON: MISSED No carboxylic acid '
               'group found\n'
               ' * SMILES: '
               'O(C(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C[C@H](O)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC '
               'NAME: DG(20:3n9/0:0/18:3n3) REASON: MISSED Carbon count (41) '
               'outside typical range for icosanoids\n'
               ' * SMILES: '
               'O(C(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C[C@@H](O)COC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC '
               'NAME: DG(18:3n3/0:0/20:4n3) REASON: MISSED Carbon count (41) '
               'outside typical range for icosanoids\n'
               ' * SMILES: '
               'OC1C(C(C(O)C1)CCCCCCC(OC(C[N+](C)(C)C)CC([O-])=O)=O)CCC(O)CCCCC '
               'NAME: '
               '7-[(1R,2R,3R,5S)-3,5-Dihydroxy-2-[(3S)-3-hydroxyoctyl]cyclopentyl]heptanoylcarnitine '
               'REASON: MISSED Carbon count (27) outside typical range for '
               'icosanoids\n'
               ' * SMILES: '
               'O(C(=O)CCC/C=C\\CCCCCCCC)C[C@@H](O)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC '
               'NAME: DG(14:1n5/0:0/20:5n3) REASON: MISSED Carbon count (37) '
               'outside typical range for icosanoids\n'
               ' * SMILES: '
               'O/1[C@@]2([C@@]([C@H]([C@H](O)C2)/C=C/[C@@H](O)CCCCC)(C\\C1=C\\CCCC(=O)NCCO)[H])[H] '
               'NAME: PGI2-EA REASON: MISSED No carboxylic acid group found\n'
               ' * SMILES: '
               'O[C@@H]1[C@@H]([C@H](C(=O)C1)/C=C/[C@@H](O)CCCCC)C/C=C\\CCCC(OC)=O '
               'NAME: Prostaglandin D2 methyl ester REASON: MISSED No '
               'carboxylic acid group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)NC4=CC=C(C=C4)OC)C5=CC=CC=C5N2C)[C@H](C)CO',
                                     'name': 'LSM-30659',
                                     'reason': 'No characteristic icosanoid '
                                               'carbon chain found'},
                                 {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'Too few oxygen atoms for an '
                                               'icosanoid'},
                                 {   'smiles': 'ClC=1C(O)=C2O[C@@]([C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC)[C@@H](O)C(=C2)C1)C(C)=C)C(=O)N3[C@@H](C=CC3)C(=O)N/C(=C(/CC)\\C)/C(=O)N/C(=C/C(O)=O)/C(O)=O)(CC)C',
                                     'name': 'Phomopsin A',
                                     'reason': 'No characteristic icosanoid '
                                               'carbon chain found'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'No characteristic icosanoid '
                                               'carbon chain found'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'No characteristic icosanoid '
                                               'carbon chain found'},
                                 {   'smiles': 'O(C1C(O)C(OC(OC=2C3=C(C=CC2)C=C(C(=C3O)C(=O)C)C)C1O)CO)C4OC(C(O)C(O)C4O)CO',
                                     'name': 'Orientaloside',
                                     'reason': 'No characteristic icosanoid '
                                               'carbon chain found'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'No characteristic icosanoid '
                                               'carbon chain found'},
                                 {   'smiles': 'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=NC=C3)O[C@@H]1CN(C)C(=O)NC4=CC5=C(C=C4)OCO5)[C@@H](C)CO',
                                     'name': '3-(1,3-benzodioxol-5-yl)-1-[[(4S,5S)-2-[(2S)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-8-(2-pyridin-4-ylethynyl)-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-1-methylurea',
                                     'reason': 'No characteristic icosanoid '
                                               'carbon chain found'},
                                 {   'smiles': '[H][C@@]12C[C@]1(C)[C@H](O)C[C@]1(C)CC[C@@]3(C)CC[C@H](C(C)C)[C@]3([H])[C@@]21[H]',
                                     'name': '13-epi-neoverrucosan-5beta-ol',
                                     'reason': 'Too few oxygen atoms for an '
                                               'icosanoid'},
                                 {   'smiles': 'O([C@@H]1[C@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]3CO)O)O[C@@H]([C@H]1O)CO[C@H]4O[C@@H]([C@@H](O)[C@H](O)[C@@H]4O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9O)CO[C@]%10(O[C@H]([C@H](NC(=O)CO)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO',
                                     'name': '(2R,4S,5R,6R)-2-[[(2R,3R,4S,5R,6S)-6-[(2R,3S,4R,5R,6S)-5-Acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2S,3S,4S,5R,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-6-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-4-hydroxy-6-(hydroxymethyl)-5-[(2S,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4-hydroxy-5-[(2-hydroxyacetyl)amino]-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No characteristic icosanoid '
                                               'carbon chain found'}],
    'sample_false_negatives': [   {   'smiles': 'O=C1OC[C@]23[C@H]4C(=O)C=5C=6C(C(=O)N(C5)[C@H](C(=O)O)CC(C)C)=C(O)C=C(C6[C@@]4(OC)[C@H]([C@@H]2OC(=O)C)C(C7=C3C1=C(O)C=C7C)=O)C',
                                      'name': 'Talauxin L',
                                      'reason': 'No characteristic icosanoid '
                                                'carbon chain found'},
                                  {   'smiles': 'O=C1OC[C@]23[C@H]4C(=O)C=5C=6C(C(=O)N(C5)[C@H](C(=O)O)C(C)C)=C(O)C=C(C6[C@@]4(OC)[C@H]([C@@H]2OC(=O)C)C(C7=C3C1=C(O)C=C7C)=O)C',
                                      'name': 'Talauxin V',
                                      'reason': 'No characteristic icosanoid '
                                                'carbon chain found'},
                                  {   'smiles': 'O=C1OC[C@]23[C@H]4C(=O)C=5C=6C(C(=O)N(C5)[C@H](C(=O)O)[C@H](CC)C)=C(O)C=C(C6[C@@]4(OC)[C@H]([C@@H]2OC(=O)C)C(C7=C3C1=C(O)C=C7C)=O)C',
                                      'name': 'Talauxin I',
                                      'reason': 'No characteristic icosanoid '
                                                'carbon chain found'},
                                  {   'smiles': 'O=C1OC[C@]23[C@H]4C(=O)C=5C=6C(C(=O)N(C5)[C@H](C(=O)O)CCC(=O)O)=C(O)C=C(C6[C@@]4(OC)[C@H]([C@@H]2OC(=O)C)C(C7=C3C1=C(O)C=C7C)=O)C',
                                      'name': 'Talauxin E',
                                      'reason': 'No characteristic icosanoid '
                                                'carbon chain found'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the issues with the previous code and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Carbon Count Range: The previous code was too restrictive '
                 '(15-25 carbons). Looking at the false negatives, many valid '
                 'icosanoids are actually part of larger molecules like '
                 'diglycerides (DG) and phosphatidylserines (PS). We should '
                 'either:\n'
                 '   - Ignore the carbon count restriction entirely, or\n'
                 '   - Look for a C20 backbone specifically within larger '
                 'molecules\n'
                 '\n'
                 '2. Carboxylic Acid Requirement: Some valid icosanoids (like '
                 "PGI2-EA and prostaglandin esters) don't have free carboxylic "
                 "acid groups because they're derivatized. We should make this "
                 'requirement optional.\n'
                 '\n'
                 '3. Missing Pattern Recognition: We should add more specific '
                 'patterns for:\n'
                 '   - Prostaglandin-like structures\n'
                 '   - Characteristic spacing of double bonds\n'
                 '   - Common modifications like hydroperoxy groups (-OOH)\n'
                 '\n'
                 '4. Molecular Weight: The current range (250-500) is too '
                 'restrictive for modified icosanoids.\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 709,
    'num_false_positives': 100,
    'num_true_negatives': 292,
    'num_false_negatives': 4,
    'num_negatives': None,
    'precision': 0.8763906056860321,
    'recall': 0.9943899018232819,
    'f1': 0.9316688567674112,
    'accuracy': 0.9058823529411765,
    'negative_predictive_value': 0.9864864864864865}