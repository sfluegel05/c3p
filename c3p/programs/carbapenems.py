"""
Classifies: CHEBI:46633 carbapenems
"""
from rdkit import Chem

def is_carbapenems(smiles: str):
    """
    Determines if a molecule is a carbapenem based on its SMILES string.
    A carbapenem contains a bicyclic ring structure with a β-lactam 
    fused to a five-membered ring and might often include a sulfur atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a carbapenem, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for carbapenem core structures
    # 4-membered ring (beta-lactam) fused to a 5-membered ring with typical substitutions
    carbapenem_patterns = [
        Chem.MolFromSmarts("C1C(=O)N2CC(C)C12"),  # Simplified for illustration
        Chem.MolFromSmarts("C1[C@H]2N(C1=O)C=C[C@H]2C")  # Bicyclic core
    ]

    matches_pattern = any(mol.HasSubstructMatch(pattern) for pattern in carbapenem_patterns)
    
    if not matches_pattern:
        return False, "Does not match any carbapenem structural core motifs"

    # Check for the presence of sulfur which is common in many carbapenems
    has_sulfur = any(atom.GetAtomicNum() == 16 for atom in mol.GetAtoms())

    if not has_sulfur:
        return True, "Matches carbapenem core structure but lacks typical sulfur atom"

    return True, "Matches carbapenem core structure with possible sulfur atom"

# Ensure the SMARTS patterns and additional checks comprehensively capture the desired chemical class.


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:46633',
                          'name': 'carbapenems',
                          'definition': 'The class of beta-lactam antibiotics '
                                        'that whose members have a carbapenem '
                                        'skeleton which is variously '
                                        'substituted at positions 3, 4, and 6.',
                          'parents': [   'CHEBI:27171',
                                         'CHEBI:27933',
                                         'CHEBI:38101'],
                          'xrefs': [   'PMID:18076336',
                                       'PMID:19199922',
                                       'PMID:19610642',
                                       'PMID:19617229',
                                       'PMID:19929840',
                                       'PMID:21859938',
                                       'Wikipedia:Carbapenem'],
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
               '[H][C@]12[C@@H](C)C(SC3Cn4cnc[n+]4C3)=C(N1C(=O)[C@]2([H])[C@@H](C)O)C([O-])=O '
               'NAME: biapenem REASON: MISSED Does not match the carbapenem '
               'structural core motif\n'
               ' * SMILES: '
               'S(C1[C@@H](N2C(=O)C(C2C1)[C@@H](O)C)C(=O)O)CCNC(=O)C NAME: '
               '17927 D REASON: MISSED Does not match the carbapenem '
               'structural core motif\n'
               ' * SMILES: '
               'S(C1=C(N2C(=O)C(C2C1C)C(O)C)C(=O)O)C3CNC(C3)C(=O)N(C)C.[Na+].[Na+].O=C([O-])[O-] '
               'NAME: Meropenem sodium carbonate REASON: MISSED Does not match '
               'the carbapenem structural core motif\n'
               ' * SMILES: '
               'S(C=1C[C@]2(N(C(=O)[C@]2([C@@H](OS(O)(=O)=O)C)[H])C1C(O)=O)[H])CCNC(=O)CCNC(=O)[C@H](O)C(CO)(C)C '
               'NAME: OA-6129 C REASON: MISSED Does not match the carbapenem '
               'structural core motif\n'
               ' * SMILES: '
               'S(=O)(C1=C(N2C(=O)[C@@H]([C@H]2C1)C(O)(C)C)C(=O)O)CCNC(=O)C '
               'NAME: Carpetimycin C REASON: MISSED Does not match the '
               'carbapenem structural core motif\n'
               ' * SMILES: '
               'C(OCOC(C(C)(C)C)=O)(=O)C=1N2[C@]([C@H](C1SC3CN(C3)C4=NCCS4)C)([C@](C2=O)([C@@H](C)O)[H])[H] '
               'NAME: tebipenem pivoxil REASON: MISSED Does not match the '
               'carbapenem structural core motif\n'
               ' * SMILES: '
               'N12[C@@]([C@]([C@@H](OS(O)(=O)=O)C)(C1=O)[H])(CC(=C2C(O)=O)SCCN)[H] '
               'NAME: Deacetylepithienamycin F REASON: MISSED Does not match '
               'the carbapenem structural core motif\n'
               ' * SMILES: '
               'O.[C@@H]1(SC2=C(N3C([C@]([C@]3([C@H]2C)[H])([H])[C@@H](C)O)=O)C(O)=O)CN[C@@H](C1)CNS(N)(=O)=O '
               'NAME: Doripenem hydrate REASON: MISSED Does not match the '
               'carbapenem structural core motif\n'
               ' * SMILES: '
               'S(C1=C(N2C(=O)\\C(\\[C@@H]2C1)=C(\\CO)/C)C(=O)O)/C=C/NC(=O)C '
               'NAME: Asparenomycin C REASON: MISSED Does not match the '
               'carbapenem structural core motif\n'
               ' * SMILES: [H][C@]12CC=C(N1C(=O)C2)C(O)=O NAME: '
               '1-carbapenem-3-carboxylic acid REASON: MISSED Does not match '
               'the carbapenem structural core motif\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'Does not match any carbapenem '
                                               'structural core motifs'}],
    'sample_false_negatives': [   {   'smiles': 'S(C1[C@@H](N2C(=O)C(C2C1)[C@@H](O)C)C(=O)O)CCNC(=O)C',
                                      'name': '17927 D',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': 'S(C=1C[C@]2(N(C(=O)[C@]2([C@@H](OS(O)(=O)=O)C)[H])C1C(O)=O)[H])CCNC(=O)CCNC(=O)[C@H](O)C(CO)(C)C',
                                      'name': 'OA-6129 C',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': 'S(=O)(C1=C(N2C(=O)[C@@H]([C@H]2C1)C(O)(C)C)C(=O)O)CCNC(=O)C',
                                      'name': 'Carpetimycin C',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': 'N12[C@@]([C@]([C@@H](OS(O)(=O)=O)C)(C1=O)[H])(CC(=C2C(O)=O)SCCN)[H]',
                                      'name': 'Deacetylepithienamycin F',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': 'S(C1=C(N2C(=O)\\C(\\[C@@H]2C1)=C(\\CO)/C)C(=O)O)/C=C/NC(=O)C',
                                      'name': 'Asparenomycin C',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': '[H][C@]12CC=C(N1C(=O)C2)C(O)=O',
                                      'name': '1-carbapenem-3-carboxylic acid',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': 'C[C@H]([C@H]1[C@H]2CC(=C(N2C1=O)C(=O)O)SCCN=CN)O',
                                      'name': '(5R,6R)-3-[2-(aminomethylideneamino)ethylthio]-6-[(1R)-1-hydroxyethyl]-7-oxo-1-azabicyclo[3.2.0]hept-2-ene-2-carboxylic '
                                              'acid',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': 'S(C1=C(N2C(=O)[C@@H]([C@H]2C1)CC)C(=O)O)CCNC(=O)CCNC(=O)[C@H](O)C(CO)(C)C',
                                      'name': 'OA-6129 A',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': 'S(=O)(=O)(OC(C1C(=O)N2C(C(=O)O)=C(CC12)S(=O)/C=C/NC(=O)C)C)O',
                                      'name': 'MM 4550',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'},
                                  {   'smiles': 'S(=O)(C1=C(N2C(=O)[C@H]([C@H]2C1)C(O)(C)C)C(=O)O)/C=C/NC(=O)C',
                                      'name': 'Carpetimycin A',
                                      'reason': 'Does not match any carbapenem '
                                                'structural core motifs'}],
    'attempt': 3,
    'reasoning': "To improve the classification of carbapenems, let's first "
                 'analyze the reasons for the failure in the previous code and '
                 'address them:\n'
                 '\n'
                 '1. **SMARTS Pattern Complexity**:\n'
                 '   The SMARTS pattern used was overly simplistic and did not '
                 'accurately capture the specific characteristics of '
                 'carbapenems. Carbapenems have distinct features, especially '
                 'the β-lactam fused with a five-membered ring, along with '
                 'other possible substituents like sulfur.\n'
                 '\n'
                 '2. **Substituents and Variability**:\n'
                 '   Carbapenems may have various substituents attached to the '
                 'bicyclic core, such as sulfur or hydroxyl groups, which the '
                 'previous pattern did not account for. It’s crucial to '
                 'capture these variations.\n'
                 '\n'
                 '3. **Additional Structural Motifs**:\n'
                 '   Carbapenems might have other features like the '
                 'orientation of additional groups, which should be captured '
                 'in the pattern.\n'
                 '\n'
                 '### Improved Approach\n'
                 '\n'
                 '- **Develop a more precise SMARTS pattern**:\n'
                 '  Identify SMARTS patterns that reflect the core carbapenem '
                 'structure — the bicyclic core structure of a '
                 '1β-methylcarbapenem with β-lactam and five-membered rings. '
                 'For example, starting with known structural motifs that '
                 'identify the common parts of carbapenems.\n'
                 '\n'
                 '- **Multiple Patterns**:\n'
                 '  Utilize multiple patterns to cover the structural '
                 'diversity within carbapenems.\n'
                 '\n'
                 '- **Check Additional Substituents**:\n'
                 '  Verify the presence of sulfur, hydroxyl groups, or other '
                 'characteristic functional groups that are frequently found '
                 'in carbapenems.\n'
                 '\n'
                 "Here's the revised code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 12,
    'num_false_positives': 2,
    'num_true_negatives': 142252,
    'num_false_negatives': 34,
    'num_negatives': None,
    'precision': 0.8571428571428571,
    'recall': 0.2608695652173913,
    'f1': 0.4,
    'accuracy': 0.999747013352073,
    'negative_predictive_value': 0.9997610446565368}