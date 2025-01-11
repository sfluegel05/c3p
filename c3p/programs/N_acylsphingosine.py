"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
"""
Classifies: CHEBI:37548 N-acylsphingosine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    N-acylsphingosines are composed of sphingosine having a fatty acyl group 
    attached to the nitrogen.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acylsphingosine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for key structural features:
    
    # 1. Amide group (N-C(=O)-) - more general pattern
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"
    
    # 2. Core sphingosine structure with flexibility for stereochemistry
    # [CH2OH]-[CH]-[NH]-[C]=O backbone with nearby OH and double bond
    sphingosine_core = Chem.MolFromSmarts("[CH2][OH].[CH]([NH])[CH]([OH])*.*=*")
    if not mol.HasSubstructMatch(sphingosine_core):
        return False, "Missing core sphingosine structure"

    # 3. Check for required functional groups
    
    # Primary alcohol (CH2-OH)
    primary_alcohol = Chem.MolFromSmarts("[CH2][OH]")
    if not mol.HasSubstructMatch(primary_alcohol):
        return False, "Missing primary alcohol group"
    
    # Secondary alcohol
    secondary_alcohol = Chem.MolFromSmarts("[CH]([#6])[OH]")
    if not mol.HasSubstructMatch(secondary_alcohol):
        return False, "Missing secondary alcohol group"
    
    # Double bond in chain
    alkene = Chem.MolFromSmarts("C=C")
    if not mol.HasSubstructMatch(alkene):
        return False, "Missing double bond"

    # Count carbons and check molecular weight
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 16:  # Minimum for basic N-acylsphingosine
        return False, "Carbon count too low for N-acylsphingosine"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250:
        return False, "Molecular weight too low for N-acylsphingosine"

    # Count nitrogens - should have exactly one (in the amide)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count != 1:
        return False, f"Should have exactly 1 nitrogen, found {n_count}"

    # Count oxygens - should have at least 3 (2 OH groups + 1 C=O)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 3:
        return False, "Insufficient oxygen count for N-acylsphingosine"

    # Verify long chain nature
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Too few rotatable bonds for N-acylsphingosine structure"

    return True, "Contains sphingosine backbone with N-acyl group and characteristic structural features"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:52639',
                          'name': 'N-acylsphingosine',
                          'definition': 'The parent compounds of the ceramide '
                                        'family, composed of sphingosine '
                                        'having an unspecified fatty acyl '
                                        'group attached to the nitrogen.',
                          'parents': ['CHEBI:82830', 'CHEBI:83273'],
                          'xrefs': [   'KEGG:C00195',
                                       'LIPID_MAPS_class:LMSP0201',
                                       'PMID:1378088',
                                       'PMID:7542630'],
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
               'CCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: N-2-hydroxystearoylsphingosine REASON: MISSED Missing '
               'characteristic sphingosine backbone with double bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)CCCCCNc1ccc([N+]([O-])=O)c2nonc12 '
               'NAME: '
               'N-{6-[(7-nitro-2,1,3-benzoxadiazol-4-yl)amino]hexanoyl}sphingosine '
               'REASON: MISSED Missing characteristic sphingosine backbone '
               'with double bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: N-2-hydroxylignoceroylsphingosine REASON: MISSED Missing '
               'characteristic sphingosine backbone with double bond\n'
               ' * SMILES: '
               '[C@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCOC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)CO '
               'NAME: N-(omega-linoleoyloxydotricontanoyl)sphingosine REASON: '
               'MISSED Missing characteristic sphingosine backbone with double '
               'bond\n'
               ' * SMILES: '
               '[C@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCOC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)CO '
               'NAME: N-(omega-linoleoyloxy-octacosanoyl)sphingosine REASON: '
               'MISSED Missing characteristic sphingosine backbone with double '
               'bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)C(O)CCCC NAME: '
               'N-(2-hydroxyhexanoyl)sphingosine REASON: MISSED Missing '
               'characteristic sphingosine backbone with double bond\n'
               ' * SMILES: CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(C)=O NAME: '
               'N-acetylsphingosine REASON: MISSED Missing characteristic '
               'sphingosine backbone with double bond\n'
               ' * SMILES: CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)CCCCC '
               'NAME: N-(hexanoyl)sphing-4-enine REASON: MISSED Missing '
               'characteristic sphingosine backbone with double bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCCCC NAME: '
               'N-dodecanoylsphingosine REASON: MISSED Missing characteristic '
               'sphingosine backbone with double bond\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCCCCC '
               'NAME: N-(2-hydroxyhexacosanoyl)sphingosine REASON: MISSED '
               'Missing characteristic sphingosine backbone with double bond\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'Missing double bond'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Missing double bond'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'Missing core sphingosine '
                                               'structure'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'Missing core sphingosine '
                                               'structure'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No amide group found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)CCCCCNc1ccc([N+]([O-])=O)c2nonc12',
                                      'name': 'N-{6-[(7-nitro-2,1,3-benzoxadiazol-4-yl)amino]hexanoyl}sphingosine',
                                      'reason': 'Should have exactly 1 '
                                                'nitrogen, found 5'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](COP([O-])(=O)OCC[N+](C)(C)C)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                      'name': 'N-(2-hydroxyhexacosanoyl)sphingosine',
                                      'reason': 'Missing core sphingosine '
                                                'structure'}],
    'attempt': 1,
    'reasoning': 'Let me analyze what went wrong and propose improvements:\n'
                 '\n'
                 '1. The sphingosine pattern was too specific:\n'
                 '- The previous pattern '
                 '`[CX4H2][OX2H][CX4H]([NX3H])[CX4H]([OX2H])[CX3H]=[CX3H]` was '
                 'too restrictive\n'
                 '- It required exact hydrogen counts and specific '
                 'connectivity that might not match all valid structures\n'
                 '- Some of the examples show variations in stereochemistry '
                 'and substitution patterns\n'
                 '\n'
                 '2. The alkene pattern was also too specific:\n'
                 '- The pattern '
                 '`[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX2]=[CX2]~[CX4,CX3]` was '
                 'looking for a specific chain length\n'
                 '- N-acylsphingosines can have varying chain lengths and '
                 'double bond positions\n'
                 '\n'
                 '3. The amide pattern needs to be more flexible:\n'
                 '- Should account for substituted nitrogens\n'
                 '- Should allow for different types of acyl groups\n'
                 '\n'
                 "Here's my improved version that should better handle the "
                 'structural variations:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 25,
    'num_false_positives': 100,
    'num_true_negatives': 43487,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.2,
    'recall': 0.9259259259259259,
    'f1': 0.32894736842105265,
    'accuracy': 0.9976613014169762,
    'negative_predictive_value': 0.9999540113591943}