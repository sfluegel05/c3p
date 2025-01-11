"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine contains a nitrogen atom bonded to exactly one carbon group R and
    two hydrogen atoms, fitting the formula: RNH2.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for primary amine, considering implicit hydrogens
    primary_amine_pattern = Chem.MolFromSmarts("[NX3;H2;!$(NC=O)]")

    # Check if structure matches the primary amine SMARTS pattern
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains a primary amine group (RNH2 structure identified)"

    return False, "No primary amine group found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:32877',
                          'name': 'primary amine',
                          'definition': 'A compound formally derived from '
                                        'ammonia by replacing one hydrogen '
                                        'atom by a hydrocarbyl group.',
                          'parents': ['CHEBI:32952', 'CHEBI:50994'],
                          'xrefs': [   'KEGG:C00375',
                                       'KEGG:C00893',
                                       'KEGG:C02580'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: O(C=1[N+](=NC=C(C1)N)C2=CC=CC=C2)C '
               'NAME: amezinium REASON: MISSED No primary amine group found\n'
               ' * SMILES: [O-]/[N+](/N(CCCN)CCC)=N\\O NAME: PAPA NONOate '
               'REASON: MISSED No primary amine group found\n'
               ' * SMILES: CO\\C=C\\C(=O)C1=NC2=C3C(C=NC3=C(N)C=C2NC(C)=O)=C1 '
               'NAME: lymphostin REASON: MISSED No primary amine group found\n'
               ' * SMILES: C[C@]12CC3CC(C1)(C[C@@](C3)(C2)C)N NAME: '
               '(3R,5S)-3,5-dimethyl-1-adamantanamine REASON: MISSED No '
               'primary amine group found\n'
               ' * SMILES: CC(C)(C)NC[C@H](O)c1cc(Cl)c(N)c(Cl)c1 NAME: '
               '(R)-clenbuterol REASON: MISSED No primary amine group found\n'
               ' * SMILES: Nc1ccccc1 NAME: aniline REASON: MISSED No primary '
               'amine group found\n'
               ' * SMILES: C1=CC=C(C(=C1)CCN)S(=O)(=O)F NAME: '
               '2-(2-aminoethyl)benzenesulfonyl fluoride REASON: MISSED No '
               'primary amine group found\n'
               ' * SMILES: CCCCCCCCN NAME: octan-1-amine REASON: MISSED No '
               'primary amine group found\n'
               ' * SMILES: C1=CC(=C(C=C1O)O)CCN NAME: '
               '4-(2-aminoethyl)benzene-1,3-diol REASON: MISSED No primary '
               'amine group found\n'
               ' * SMILES: NCCC=C NAME: 3-buten-1-amine REASON: MISSED No '
               'primary amine group found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'P(OCC(OC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCC/C=C\\C/C=C\\C/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O',
                                     'name': 'Pe-nme2(20:3(5Z,8Z,11Z)/22:6(4Z,7Z,10Z,13Z,16Z,19Z))',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'O1[C@]2([C@@H](O)[C@H](O)[C@@H](O)[C@@]1(OCC[C@H](CCC=C(C(OC[C@]3(O)[C@@H](O)[C@@](OC3)(OC2)[H])=O)C)C)[H])[H]',
                                     'name': 'Urceolide',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@@H]2[C@@H](CO)O[C@H](O[C@@H]3[C@@H](CO)O[C@H](O[C@@H]4[C@@H](CO)O[C@H](O)[C@H](O)[C@H]4O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@@H]1O',
                                     'name': 'alpha-maltotetraose',
                                     'reason': 'No primary amine group found'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'No primary amine group found'}],
    'sample_false_negatives': [   {   'smiles': 'CC=CC1=CC=C2[C@@H]3[C@@H]([C@H]([C@@H](N3CCC4=CC=CC=C4)CN2C1=O)CO)C(=O)N(C)C',
                                      'name': 'LSM-14647',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'ClC1=C(Cl)C=CC(=C1)CCN(C[C@@H]2N(CCCC2)C)C',
                                      'name': '2-Azabutane, '
                                              '1-[1-methylazacyclohex-2-yl]-2-methyl-4-[3,4-dichlorophenyl]-',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'N(CCC1=CC=CC=C1)(C)C',
                                      'name': 'N,N-Dimethylphenethylamine',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'S=C(NC1CCCCC1)N1CCC(CC1)c1c[nH]cn1',
                                      'name': 'Thioperamide',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'ClC=1C(CCN[C@@H]2[C@]3(O[C@@]4(O[C@H](OC)[C@@]3([C@@](C2)([C@@H]4C)[H])[H])[H])C)=CC=CC1',
                                      'name': '4-(2-(2-Chlorophenyl)ethyl)amino-8-methoxy-3,10-dimethyl-2,9-dioxatricyclo(4,3,1,0(3,7))decane',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'CN1CCN(CC1)CCC2=CC(=C(C=C2)Cl)Cl',
                                      'name': '1-[2-(3,4-dichlorophenyl)ethyl]-4-methylpiperazine',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'CN[C@@H]1[C@H](O)[C@H](NC)[C@H]2O[C@]3(O)[C@@H](O[C@H](C)C=C3OC)O[C@@H]2[C@H]1O',
                                      'name': 'Spenolimycin',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'CC(C)(C)NC(=O)COC1=C(C=C(C=C1)CNCCC2=CC=CC=C2)OC',
                                      'name': 'N-tert-butyl-2-[2-methoxy-4-[(2-phenylethylamino)methyl]phenoxy]acetamide',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'COC(=O)[C@H]1[C@@H]([C@H]2CN3C(=O)C=CC=C3[C@@H]1N2CCC4=CC=CC=C4)CO',
                                      'name': 'LSM-13982',
                                      'reason': 'No primary amine group found'},
                                  {   'smiles': 'COC(=O)[C@@H]1[C@H]([C@@H]2CN3C(=O)C=CC=C3[C@H]1N2CCC4=CC=CC=C4)CO',
                                      'name': 'LSM-12804',
                                      'reason': 'No primary amine group '
                                                'found'}],
    'attempt': 1,
    'reasoning': 'To classify a primary amine correctly, we need to ensure our '
                 'pattern recognition specifically identifies structures where '
                 'the nitrogen atom is bonded to exactly one carbon atom '
                 '(alkyl or aryl) and two hydrogen atoms, fitting the formula '
                 'RNH2. The original code attempted to do this by analyzing '
                 "nitrogen atoms' neighbors, but it might have failed to "
                 'account for various bonding contexts or configurations '
                 'represented in SMILES strings.\n'
                 '\n'
                 'Here are potential issues and improvements:\n'
                 "1. **Hydrogen atoms in SMILES strings**: SMILES doesn't "
                 'always explicitly list hydrogen atoms unless they are part '
                 'of functional groups. In many cases, implicit hydrogens are '
                 'assumed based on atom valences.\n'
                 '2. **Ring systems and different attachment points**: Our '
                 'previous pattern might not account for nitrogen atoms within '
                 'ring systems or attached via complex structures.\n'
                 '3. **Incorrect pattern logic**: SMILES with implicit details '
                 '(like aromatic systems) can affect bondage checks.\n'
                 '\n'
                 'Revise the logic to better capture primary amines by using '
                 'SMARTS pattern matching, which accounts for implicit '
                 "hydrogens as per the molecule's valency rules, thereby "
                 'improving the recognition of primary amines in various '
                 'structural contexts.\n'
                 '\n'
                 "Here's the revised code with SMARTS pattern matching:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 101,
    'num_false_positives': 100,
    'num_true_negatives': 616,
    'num_false_negatives': 18,
    'num_negatives': None,
    'precision': 0.5024875621890548,
    'recall': 0.8487394957983193,
    'f1': 0.63125,
    'accuracy': 0.858682634730539,
    'negative_predictive_value': 0.9716088328075709}