"""
Classifies: CHEBI:20060 3-hydroxy fatty acyl-CoA
"""
from rdkit import Chem

def is_3_hydroxy_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acyl-CoA based on its SMILES string.
    A 3-hydroxy fatty acyl-CoA has a 3-hydroxy group on the fatty acid and a CoA moiety.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if a 3-hydroxy fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a 3-hydroxy group attached to a fatty acyl chain
    hydroxy_pattern = Chem.MolFromSmarts("[C](O)CC(=O)")
    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "No 3-hydroxy group on fatty acyl chain found"
    
    # Simplified Coenzyme A detection focusing on key elements
    # e.g., include [SC](=O)C(=O)NCC motif as part of CoA
    coa_pattern = Chem.MolFromSmarts("C(=O)NCCSC(=O)")
    if not mol.HasSubstructMatch(coa_pattern):
        return False, "No Coenzyme A moiety found"

    return True, "Identified as 3-hydroxy fatty acyl-CoA with CoA moiety"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:20060',
                          'name': '3-hydroxy fatty acyl-CoA',
                          'definition': 'A hydroxy fatty acyl-CoA that results '
                                        'from the formal condensation of the '
                                        'thiol group of coenzyme A with the '
                                        'carboxy group of any 3-hydroxy fatty '
                                        'acid.',
                          'parents': ['CHEBI:61902'],
                          'xrefs': [   'PMID:12106015',
                                       'PMID:1778900',
                                       'PMID:20583174',
                                       'PMID:20670938',
                                       'PMID:20923481',
                                       'PMID:21502722',
                                       'PMID:7552767'],
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
               'CCCCCCCCCCC[C@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (S)-3-hydroxytetradecanoyl-CoA REASON: MISSED No '
               'appropriate 3-hydroxy group in fatty acid chain found\n'
               ' * SMILES: '
               'CCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-hydroxydecanoyl-CoA REASON: MISSED No appropriate '
               '3-hydroxy group in fatty acid chain found\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: '
               '(3R,10Z,13Z,16Z,19Z,22Z,25Z)-3-hydroxyoctacosahexaenoyl-CoA '
               'REASON: MISSED No appropriate 3-hydroxy group in fatty acid '
               'chain found\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (3R,11Z,14Z,17Z,20Z)-3-hydroxyhexacosatetraenoyl-CoA '
               'REASON: MISSED No appropriate 3-hydroxy group in fatty acid '
               'chain found\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC[C@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: '
               '(3S,6Z,9Z,12Z,15Z,18Z,21Z)-3-hydroxytetracosahexaenoyl-CoA '
               'REASON: MISSED No appropriate 3-hydroxy group in fatty acid '
               'chain found\n'
               ' * SMILES: '
               '[C@@H]1(N2C3=C(C(=NC=N3)N)N=C2)O[C@H](COP(OP(OCC(C)([C@H](C(NCCC(NCCSC(=O)C[C@@H](O)CCCCCCCCCCCCCCCCC)=O)=O)O)C)(=O)O)(=O)O)[C@H]([C@H]1O)OP(O)(O)=O '
               'NAME: (3S)-3-hydroxyicosanoyl-CoA REASON: MISSED No '
               'appropriate 3-hydroxy group in fatty acid chain found\n'
               ' * SMILES: '
               'C[C@](O)(CC(O)=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (3S)-3-hydroxy-3-methylglutaryl-CoA REASON: MISSED No '
               'appropriate 3-hydroxy group in fatty acid chain found\n'
               ' * SMILES: '
               'CC(C)(COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)[C@@H](O)C(=O)NCCC(=O)NCCSC(=O)C1=CCCCC1O '
               'NAME: 6-hydroxycyclohex-1-ene-1-carbonyl-CoA REASON: MISSED No '
               'appropriate 3-hydroxy group in fatty acid chain found\n'
               ' * SMILES: '
               'CC(C)(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-hydroxyisovaleryl-CoA REASON: MISSED No appropriate '
               '3-hydroxy group in fatty acid chain found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCC(O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-hydroxytetracosanoyl-CoA REASON: MISSED No appropriate '
               '3-hydroxy group in fatty acid chain found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No 3-hydroxy group on fatty '
                                               'acyl chain found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No 3-hydroxy group on fatty '
                                               'acyl chain found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No 3-hydroxy group on fatty '
                                               'acyl chain found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No 3-hydroxy group on fatty '
                                               'acyl chain found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No 3-hydroxy group on fatty '
                                               'acyl chain found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No 3-hydroxy group on fatty '
                                               'acyl chain found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No 3-hydroxy group on fatty '
                                               'acyl chain found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No Coenzyme A moiety found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No 3-hydroxy group on fatty '
                                               'acyl chain found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'In the previous attempt, the program failed to correctly '
                 'identify 3-hydroxy fatty acyl-CoA molecules due to overly '
                 'specific substructure patterns, which missed expected '
                 'features. Specifically, the failed classification likely '
                 'results from overly rigid SMARTS patterns for the 3-hydroxy '
                 'group and Coenzyme A moiety.\n'
                 '\n'
                 'To improve this, we can take the following steps:\n'
                 '\n'
                 '1. **Hydroxy Group Detection**: The pattern '
                 '`[C@@H](O)[CH2][CH2][CX3](=O)` may not cover all structural '
                 'variations for 3-hydroxy fatty acids. We should generalize '
                 'the pattern to account for longer and branched chains, such '
                 'as `[C](O)CC(C=O)`. Look only for the `C-O` linkage close to '
                 'carbonyl for 3-hydroxy indication.\n'
                 '\n'
                 '2. **Coenzyme A Moiety Detection**: Coenzyme A is complex '
                 'but central features, such as a pantetheine phosphate '
                 'structure, characterize it. We should broaden its SMARTS '
                 'pattern to capture the essence of CoA, focusing on '
                 'functional groups rather than the complete structure, since '
                 'the exact CoA configuration might vary.\n'
                 '\n'
                 "Here's the revised code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 75,
    'num_false_positives': 100,
    'num_true_negatives': 51984,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.42857142857142855,
    'recall': 1.0,
    'f1': 0.6,
    'accuracy': 0.9980827853294734,
    'negative_predictive_value': 1.0}