"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: acyl-CoA compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester formed between coenzyme A and a carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for thioester linkage (-C(=O)S-) - this is essential for acyl-CoA
    thioester_patterns = [
        Chem.MolFromSmarts("[CX3](=[OX1])[SX2]"),
        Chem.MolFromSmarts("C(=O)S")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in thioester_patterns):
        return False, "No thioester linkage found"

    # Check for adenine base with flexible connection point
    adenine_patterns = [
        Chem.MolFromSmarts("c1nc(N)c2ncnc2n1"),
        Chem.MolFromSmarts("n1cnc2c(N)ncnc12")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in adenine_patterns):
        return False, "No adenine moiety found"

    # Check for phosphate groups with flexible patterns
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)O"),
        Chem.MolFromSmarts("OP([O,OH])(=O)[O,OH]"),
        Chem.MolFromSmarts("OP(O)(O)=O")
    ]
    phosphate_count = 0
    for pattern in phosphate_patterns:
        if pattern:
            phosphate_count += len(mol.GetSubstructMatches(pattern))
    if phosphate_count < 2:  # Allow for some variation in phosphate representation
        return False, f"Insufficient phosphate groups found"

    # Check for pantetheine portion with flexible pattern
    pantetheine_patterns = [
        Chem.MolFromSmarts("NCCC(=O)NCCS"),
        Chem.MolFromSmarts("NCCSC(=O)"),
        Chem.MolFromSmarts("NC(=O)CCNC(=O)")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in pantetheine_patterns):
        return False, "No pantetheine moiety found"

    # Check for ribose sugar with very flexible pattern
    ribose_patterns = [
        Chem.MolFromSmarts("OC1C(O)C(O)C(O)C1O"),
        Chem.MolFromSmarts("OC1CCCO1"),  # Basic furanose pattern
        Chem.MolFromSmarts("C1OC(CO)C(O)C1"),
        Chem.MolFromSmarts("C1OC(COP)C(O)C1")  # Connected to phosphate
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in ribose_patterns):
        return False, "No ribose sugar found"

    # Check for characteristic CoA features with flexible patterns
    coa_patterns = [
        Chem.MolFromSmarts("CC(C)(COP)"),  # Geminal dimethyl
        Chem.MolFromSmarts("SCCNC(=O)"),   # Core connection
        Chem.MolFromSmarts("COP(O)OP")     # Phosphate linkage
    ]
    if not all(mol.HasSubstructMatch(pattern) for pattern in coa_patterns):
        return False, "Missing key CoA structural features"

    # If all essential elements are present, it's an acyl-CoA
    return True, "Contains complete CoA moiety with thioester linkage to acyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:17984',
                          'name': 'acyl-CoA',
                          'definition': 'A thioester that results from the '
                                        'formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any carboxylic acid.',
                          'parents': [   'CHEBI:140325',
                                         'CHEBI:20706',
                                         'CHEBI:231540',
                                         'CHEBI:51277'],
                          'xrefs': [   'KEGG:C00040',
                                       'PMID:11264983',
                                       'PMID:11524729',
                                       'PMID:16495773',
                                       'PMID:21514367',
                                       'PMID:21541677'],
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
               'S(C(CCCCCCC/C=C\\C[C@@H](CCCCCC)O)=O)CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]1O[C@@H](N2C3=C(C(=NC=N3)N)N=C2)[C@@H]([C@@H]1OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O '
               'NAME: (9Z,12R)-12-hydroxyoctadec-9-enoyl-CoA REASON: MISSED No '
               'ribose sugar found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (R)-3-hydroxypalmitoyl-CoA REASON: MISSED No ribose '
               'sugar found\n'
               ' * SMILES: '
               'C(CCC)C/C=C\\C/C=C\\CC1C(C/C=C\\CCCC(SCCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP(O)(O)=O)O)(=O)O)(=O)O)(C)C)O)=O)=O)O1 '
               'NAME: 8,9-epoxy-(5Z,11Z,14Z)-icosatrienoyl-CoA REASON: MISSED '
               'No ribose sugar found\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)c1ccccc1N '
               'NAME: anthraniloyl-CoA REASON: MISSED No ribose sugar found\n'
               ' * SMILES: '
               'CC(=C)C(CCCC(O)=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-isopropenylpimeloyl-CoA REASON: MISSED No ribose sugar '
               'found\n'
               ' * SMILES: '
               'CCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (R)-3-hydroxyhexanoyl-CoA REASON: MISSED No ribose sugar '
               'found\n'
               ' * SMILES: '
               'CC(=O)CC(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-hydroxy-5-oxohexanoyl-CoA REASON: MISSED No ribose '
               'sugar found\n'
               ' * SMILES: '
               'CCCCC\\C=C\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: trans,trans-deca-2,4-dienoyl-CoA REASON: MISSED No '
               'ribose sugar found\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: '
               '(3R,10Z,13Z,16Z,19Z,22Z,25Z)-3-hydroxyoctacosahexaenoyl-CoA '
               'REASON: MISSED No ribose sugar found\n'
               ' * SMILES: O1C(CCC1)C(=O)N2CCNCC2 NAME: '
               '1-(Tetrahydro-2-furoyl)piperazine REASON: MISSED No adenine '
               'moiety found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H](C1=CC=CC=C1)NC(=O)C[C@H]2CC[C@@H]([C@@H](O2)CO)NC(=O)CN3CCOCC3',
                                     'name': '2-[(2R,5S,6R)-6-(hydroxymethyl)-5-[[2-(4-morpholinyl)-1-oxoethyl]amino]-2-oxanyl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'C[N+](C)(C)[C@@H](Cc1c[nH]c(n1)S(=O)C[C@H](NC(=O)CC[C@H]([NH3+])C([O-])=O)C([O-])=O)C([O-])=O',
                                     'name': 'N(alpha)-(L-gamma-glutamyl)-hercynyl-L-cysteine '
                                             'sulfoxide(1-)',
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'O(C1=CC=2[C@]3([C@](N(CC3)C)(N(C2C=C1)C)[H])C)C(=O)N4CCC=5C(C4)=CC=CC5',
                                     'name': 'quilostigmine',
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'O[C@@H]1[C@]23[C@@]4(N(C[C@@]([C@]2(C[C@@]4([C@]56[C@]3(CC(=O)[C@](C5)(C([C@H]6O)=C)[H])[H])[H])[H])(CC1)C)CC)[H]',
                                     'name': 'Bullatine G',
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CCC(O)=O)[C@H]3O',
                                     'name': "3'-L-glutamyl-AMP",
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'O1C2(C(C3(C(C4(C(CC3OC(=O)C)C(OC(=O)CC4)(C)C)C)CC2)C)CC15C6N(C=7C5=CC=CC7)C(=O)C(N6)C)C',
                                     'name': 'Teraspiridole C_130091',
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'O(C1[C@@H](OC(=O)C)C(O[C@@H](OC2=C(OC3=C(C2=O)C(O)=CC(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)CO)=C3CC=C(C)C)C5=CC=C(OC)C=C5)[C@H]1O)C)[C@@H]6OC[C@@H](O)[C@H](OC(=O)C)C6O',
                                     'name': 'Sempervirenoside A',
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'COC[C@]1(C(=O)C2CCN1CC2)CO',
                                     'name': '(2S)-2-(hydroxymethyl)-2-(methoxymethyl)-1-azabicyclo[2.2.2]octan-3-one',
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'Oc1c(C2CC(Cc3ccccc23)c2ccc(OCc3ccc(cc3)C(F)(F)F)cc2)c(=O)oc2ccccc12',
                                     'name': 'Flocoumafen',
                                     'reason': 'No thioester linkage found'},
                                 {   'smiles': 'O[C@H]1CC=2C(N(C=3C1=CC=CC3)C(=O)N)=CC=CC2',
                                     'name': '(S)-MHD',
                                     'reason': 'No thioester linkage found'}],
    'sample_false_negatives': [   {   'smiles': 'O1C(CCC1)C(=O)N2CCNCC2',
                                      'name': '1-(Tetrahydro-2-furoyl)piperazine',
                                      'reason': 'No thioester linkage found'},
                                  {   'smiles': 'O=NN1CCNCC1',
                                      'name': 'N-Mononitrosopiperazine',
                                      'reason': 'No thioester linkage found'}],
    'attempt': 2,
    'reasoning': 'Let me analyze what went wrong and suggest improvements:\n'
                 '\n'
                 '1. The main issue seems to be that the SMARTS patterns are '
                 'too strict, particularly for the ribose sugar. Many valid '
                 "acyl-CoAs are being missed because they don't match the "
                 'exact stereochemistry patterns.\n'
                 '\n'
                 '2. Looking at the false negatives, they all contain the core '
                 'CoA structure but were rejected mainly due to the "No ribose '
                 'sugar found" check. This suggests our ribose pattern is too '
                 'restrictive.\n'
                 '\n'
                 '3. The geminal dimethyl and CoA backbone patterns might also '
                 'be too specific, not accounting for potential variations in '
                 'representation.\n'
                 '\n'
                 "Here's my improved approach:\n"
                 '\n'
                 '1. Make the ribose pattern more flexible by:\n'
                 '   - Removing specific stereochemistry requirements\n'
                 '   - Using simpler ring patterns\n'
                 '   - Allowing for different representations of the sugar\n'
                 '\n'
                 '2. Simplify the CoA backbone check to focus on essential '
                 'connectivity\n'
                 '\n'
                 '3. Add additional validation patterns to catch different '
                 'SMILES representations\n'
                 '\n'
                 '4. Make the phosphate pattern more flexible to catch '
                 'different representations\n'
                 '\n'
                 "Here's my improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 636,
    'num_false_positives': 100,
    'num_true_negatives': 15305,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.8641304347826086,
    'recall': 0.9968652037617555,
    'f1': 0.9257641921397379,
    'accuracy': 0.9936420868914791,
    'negative_predictive_value': 0.9998693408244594}