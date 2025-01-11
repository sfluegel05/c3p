"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:37671 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is a sugar with one or more hydroxyl groups replaced by amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ring structure
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No ring structure found"
    
    # Look for at least 2 hydroxyl groups (sugar-like)
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    if len(hydroxyl_matches) < 2:
        return False, "Not enough hydroxyl groups for a sugar"

    # Check for amino groups (including N-acetyl and other substitutions)
    amino_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N=O)]")
    amino_matches = mol.GetSubstructMatches(amino_pattern)
    
    # Check for N-acetyl groups specifically
    n_acetyl_pattern = Chem.MolFromSmarts("[NX3;H0]([C&R])(C=O)")
    n_acetyl_matches = mol.GetSubstructMatches(n_acetyl_pattern)
    
    # Count total amino groups (including N-acetyl)
    total_amino = len(amino_matches) + len(n_acetyl_matches)
    if total_amino < 1:
        return False, "No amino groups found"

    # Check if amino groups are attached to ring carbons
    amino_sugar_pattern = Chem.MolFromSmarts("[C&R]-[NX3;H2,H1,H0;!$(N=O)]")
    if not mol.HasSubstructMatch(amino_sugar_pattern):
        return False, "Amino group not attached to ring carbon"

    # More flexible sugar pattern matching
    # Look for at least 3 contiguous carbons with hydroxyl or amino groups
    sugar_pattern = Chem.MolFromSmarts("[C&R][C&R][C&R]")
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar-like carbon chain found"

    # Additional check for common sugar functional groups
    sugar_functional_groups = Chem.MolFromSmarts("[C&R][OH]")
    if not mol.HasSubstructMatch(sugar_functional_groups):
        return False, "No sugar-like functional groups found"

    return True, f"Contains a sugar structure with {total_amino} amino group(s) attached"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:28963',
                          'name': 'amino sugar',
                          'definition': 'Any sugar having one or more '
                                        'alcoholic hydroxy groups replaced by '
                                        'substituted or unsubstituted amino '
                                        'groups.',
                          'parents': ['CHEBI:63299'],
                          'xrefs': [   'KEGG:C05383',
                                       'PMID:18424273',
                                       'PMID:9056391'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.038461538461538464 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CN[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O[C@H]1O[C@H]1[C@H](O[C@H]2[C@H](O)[C@@H](O)[C@H](OC(N)=O)[C@@H](O)[C@@H]2NC(N)=N)O[C@@H](C)[C@]1(O)CO '
               'NAME: Bluensomycin REASON: CORRECT Contains a sugar structure '
               'with 4 amino group(s) attached\n'
               ' * SMILES: '
               'O=C(N[C@H]1[C@H](O[C@H](CO)[C@H]([C@@H]1O)O)O[C@@]2(C=C[C@@H]([C@H](CCC=C(C)C)C)CC2)C)C '
               'NAME: Trichaspside A REASON: CORRECT Contains a sugar '
               'structure with 1 amino group(s) attached\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'O([C@H]1[C@@H](O[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O[C@@H](O[C@@H]3[C@@H](NC(=O)C)[C@@H](O[C@@H]([C@H]3O)CO)O[C@@H]4[C@@H](O)[C@@H](O[C@@H]([C@@H]4O)CO)O[C@H]5[C@H](O)[C@@H](O)[C@@H](O[C@@H]5CO)O)[C@@H]1O)CO)[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O '
               'NAME: '
               '(2S,4S,5R,6R)-5-Acetamido-2-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5R,6R)-3-acetamido-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2-[(2S,3R,4R,5S,6R)-3-acetamido-2-[(2R,3S,4S,5R,6S)-3,5-dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-5-hydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-3-hydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: MISSED No sugar-like carbon chain found\n'
               ' * SMILES: '
               'S(O[C@@H]1[C@@H](O)[C@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@H](O[C@H]3[C@H](O)[C@H](O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)[C@@H](O[C@H]7[C@H](O)[C@H](O[C@@H](O[C@H]8[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]8CO)O[C@H]9[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]9CO[C@@H]%10O[C@H]([C@@H](O)[C@@H](O)[C@@H]%10O)C)O)[C@H]7O)CO[C@H]%11O[C@@H]([C@@H](O)[C@H](O)[C@@H]%11O[C@@H]%12O[C@@H]([C@@H](O[C@@H]%13O[C@@H]([C@H](O)[C@H](O[C@]%14(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%14)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%13O)CO)[C@H](O)[C@H]%12NC(=O)C)CO)CO[C@@H]%15O[C@@H]([C@@H](O[C@@H]%16O[C@@H]([C@H](O)[C@H](O[C@]%17(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%17)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%16O)CO)[C@H](O)[C@H]%15NC(=O)C)CO)O[C@@H]3CO)O[C@@H]2CO)O[C@@H]([C@@H]1O)CO)(O)(=O)=O '
               'NAME: PubChem CID 91849816 REASON: MISSED No sugar-like carbon '
               'chain found\n'
               ' * SMILES: '
               'O1[C@@H](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)[C@H](NC(=O)C)[C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)CO)[C@H](O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)[C@H]1CO '
               'NAME: '
               'N-[(3R,4R,5S,6R)-5-[(2S,3R,4S,5S,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-6-(hydroxymethyl)-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: MISSED No sugar-like carbon chain found\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)[C@@H]4O)CO[C@H]6O[C@@H]([C@@H](O)[C@H](O)[C@@H]6O[C@@H]7O[C@@H]([C@@H](O)[C@H](O[C@@H]8O[C@H]([C@@H](O)[C@@H](O)[C@@H]8O)C)[C@H]7NC(=O)C)CO)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2S,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4S,5R,6R)-6-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-hydroxy-6-(hydroxymethyl)-4-[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-[[(2R,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxymethyl]oxan-3-yl]acetamide '
               'REASON: MISSED No sugar-like carbon chain found\n'
               ' * SMILES: '
               'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O)[C@@H]3O[C@@H]([C@@H](O[C@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O)CO)[C@H](O)[C@@H]4O)CO)[C@H](O)[C@@H]3O)CO '
               'NAME: '
               'N-[(2R,3R,4R,5S,6R)-5-[(2R,3R,4R,5S,6R)-3-Acetamido-5-[(2S,3S,4R,5S,6R)-5-[(2R,3S,4R,5S,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-[(2S,3S,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxy-3,4-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide '
               'REASON: MISSED No sugar-like carbon chain found\n'
               ' * SMILES: '
               'O([C@@H]1[C@@H](O)[C@@H](O[C@H]2[C@H](O[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)[C@H](O[C@@H](O)[C@@H]2NC(=O)C)CO)O[C@@H]([C@@H]1O)CO)[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O '
               'NAME: '
               '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3R,4R,5S,6R)-3-acetamido-2-hydroxy-6-(hydroxymethyl)-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-4-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
               'acid REASON: MISSED No sugar-like carbon chain found\n'
               ' * SMILES: '
               'CC(=O)N[C@H]1C(O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@H]5NC(C)=O)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O '
               'NAME: '
               'beta-D-GlcpNAc-(1->2)-[beta-D-GlcpNAc-(1->6)]-alpha-D-Manp-(1->6)-[beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)]-beta-D-Manp-(1->4)-beta-D-GlcpNAc-(1->4)-D-GlcpNAc '
               'REASON: MISSED No sugar-like carbon chain found\n'
               ' * SMILES: '
               'O([C@]1(O[C@H]([C@H](NC(=O)C)[C@@H](O)C1)[C@H](O)[C@H](O)CO)C(=O)NC)[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O '
               'NAME: '
               '(2S,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4S,5R,6S)-3,5-dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-4-hydroxy-N-methyl-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxamide '
               'REASON: MISSED No sugar-like carbon chain found\n'
               ' * SMILES: '
               'C[C@@H]1O[C@@H](OC[C@H]2O[C@@H](O)[C@H](NC(C)=O)[C@@H](O[C@@H]3O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]3O)[C@@H]2O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](O)[C@H](O)[C@@H]1O '
               'NAME: '
               'alpha-L-Fucp-(1->3)-[alpha-L-Fucp-(1->6)]-{beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-[beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->6)]-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)}-beta-D-GlcpNAc '
               'REASON: MISSED No sugar-like carbon chain found\n'
               ' * SMILES: O1C(C(O)C(OC(C)C(O)=O)C(NC(=O)C)C1O)CO NAME: '
               '(+)-N-Acetylmuramic acid REASON: MISSED No sugar-like carbon '
               'chain found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O1C(CCC=2C1=CC(OC3OC(C(O)C(O)C3O)C(O)=O)=C(C2OC)C=4C(=O)C=5C(OC4)=CC(O)=C(O)C5)(C)C',
                                     'name': '6-{[6-(6,7-dihydroxy-4-oxo-4H-chromen-3-yl)-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-7-yl]oxy}-3,4,5-trihydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No amino groups found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C',
                                     'name': 'Concanamycin D',
                                     'reason': 'No amino groups found'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'No amino groups found'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'No ring structure found'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'Amino group not attached to '
                                               'ring carbon'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar'},
                                 {   'smiles': '[C@@H]1(C2=CC=CC=C2)[C@@H](N)C1',
                                     'name': '(1S,2R)-tranylcypromine',
                                     'reason': 'Not enough hydroxyl groups for '
                                               'a sugar'},
                                 {   'smiles': 'OC(=O)/C=C/C#CC#CC#N',
                                     'name': '(2E)-7-Cyanohept-2-en-4,6-diynoic '
                                             'acid',
                                     'reason': 'No ring structure found'}],
    'sample_false_negatives': [   {   'smiles': 'CN[C@H]1[C@@H](C[C@H](N)[C@@H](O[C@H]2O[C@H](CN)CC[C@H]2N)[C@@H]1O)OC',
                                      'name': 'Istamycin Y0',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar'},
                                  {   'smiles': '[H]C([H])(C([H])([H])[C@]([H])(O)CN)[C@@]([H])(N)C=O',
                                      'name': 'purpurosamine C',
                                      'reason': 'No ring structure found'},
                                  {   'smiles': '[H]C(=O)[C@@H](N)[C@@H](O)[C@H](O)[C@H](O)CO',
                                      'name': 'aldehydo-D-mannosamine',
                                      'reason': 'No ring structure found'},
                                  {   'smiles': 'CNC[C@@H]1CC[C@@H](N)[C@@H](O[C@@H]2[C@@H](N)C[C@H](OC)[C@@H](NC)[C@H]2O)O1',
                                      'name': 'Sannamycin B',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar'},
                                  {   'smiles': 'CCNC[C@@H]1CC[C@@H](N)[C@@H](O[C@@H]2[C@@H](N)C[C@H](OC)[C@H]([C@H]2O)N(C)C(=O)CN)O1',
                                      'name': 'Istamycin C',
                                      'reason': 'Not enough hydroxyl groups '
                                                'for a sugar'},
                                  {   'smiles': 'N[C@H]([C@@H](O)CC(=O)C(O)=O)[C@@H](O)[C@H](O)[C@H](O)CO',
                                      'name': 'keto-neuraminic acid',
                                      'reason': 'No ring structure found'},
                                  {   'smiles': 'CC(=O)N[C@H]([C@@H](O)CC(=O)C(O)=O)[C@@H](O)[C@H](O)[C@H](O)CO',
                                      'name': 'aceneuramic acid',
                                      'reason': 'No ring structure found'},
                                  {   'smiles': '[H]C([H])(C=O)[C@@]([H])(N)[C@]([H])(O)[C@]([H])(C)O',
                                      'name': 'acosamine',
                                      'reason': 'No ring structure found'},
                                  {   'smiles': '[H][C@](C)(O)[C@@]([H])(N)[C@]([H])(O)[C@@]([H])(N)C=O',
                                      'name': 'bacillosamine',
                                      'reason': 'No ring structure found'},
                                  {   'smiles': 'C[C@H](\\C=C(C)\\C=C\\C(=O)NO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)C(=O)c1ccc(cc1)N(C)C',
                                      'name': 'trichostatin C',
                                      'reason': 'Amino group not attached to '
                                                'ring carbon'}],
    'attempt': 2,
    'reasoning': 'The previous code had several issues that led to a low F1 '
                 'score:\n'
                 '\n'
                 '1. **Overly Strict Sugar Pattern Matching**: The sugar '
                 'pattern `[C&R][C&R][C&R][C&R][C&R][C&R]` was too strict and '
                 'missed many valid amino sugars, especially those with '
                 'shorter carbon chains or different ring sizes.\n'
                 '\n'
                 '2. **Insufficient Handling of Complex Structures**: The code '
                 'failed to recognize complex amino sugars with multiple rings '
                 'or branched structures, which are common in the examples '
                 'provided.\n'
                 '\n'
                 '3. **Incorrect Handling of N-Acetyl Groups**: The code '
                 "counted N-acetyl groups separately but didn't properly "
                 'integrate them into the amino group count, leading to missed '
                 'classifications.\n'
                 '\n'
                 '4. **Lack of Flexibility in Hydroxyl Group Count**: The '
                 'requirement of at least 3 hydroxyl groups might be too '
                 'strict for some amino sugars.\n'
                 '\n'
                 "Here's an improved version of the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 2968,
    'num_false_positives': 100,
    'num_true_negatives': 1092,
    'num_false_negatives': 37,
    'num_negatives': None,
    'precision': 0.9674054758800521,
    'recall': 0.9876871880199667,
    'f1': 0.9774411328832537,
    'accuracy': 0.9673576364069574,
    'negative_predictive_value': 0.9672276350752879}