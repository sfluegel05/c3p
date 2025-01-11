"""
Classifies: CHEBI:23824 diol
"""
"""
Classifies: CHEBI:23824 diol
"""
from rdkit import Chem

def is_diol(smiles: str):
    """
    Determines if a molecule is a diol based on its SMILES string.
    A diol is a compound that contains two hydroxy groups, generally assumed to be, but not necessarily, alcoholic.
    Aliphatic diols are also called glycols.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for hydroxyl groups (OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    hydroxyl_count = len(hydroxyl_matches)

    if hydroxyl_count >= 2:
        return True, f"Molecule contains at least two hydroxyl groups"
    else:
        return False, f"Molecule contains {hydroxyl_count} hydroxyl group(s), diols must have at least two"

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23824',
                              'name': 'diol',
                              'definition': 'A compound that contains two hydroxy '
                                            'groups, generally assumed to be, but not '
                                            'necessarily, alcoholic. Aliphatic diols are '
                                            'also called glycols.',
                              'parents': []},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
        'message': None,
        'attempt': 2,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None,
        'num_true_positives': None,
        'num_false_positives': None,
        'num_true_negatives': None,
        'num_false_negatives': None,
        'num_negatives': None,
        'precision': None,
        'recall': None,
        'f1': None,
        'accuracy': None}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:23824',
                          'name': 'diol',
                          'definition': 'A compound that contains two hydroxy '
                                        'groups, generally assumed to be, but '
                                        'not necessarily, alcoholic. Aliphatic '
                                        'diols are also called glycols.',
                          'parents': ['CHEBI:26191'],
                          'xrefs': ['Wikipedia:Diol'],
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
               '[C@H]12N3C=4N=C(NC(C4N=C3[C@H]([C@]([C@H](C1)O)(O2)[H])O)=O)N '
               "NAME: 8,5'-cyclo-2'-deoxyguanosine REASON: MISSED Molecule "
               'contains 0 valid hydroxyl groups, diols must have at least '
               'two\n'
               ' * SMILES: OC(CNCC(O)C)C NAME: diisopropanolamine REASON: '
               'MISSED Molecule contains 0 valid hydroxyl groups, diols must '
               'have at least two\n'
               ' * SMILES: O[C@@H]1C=CC=C([C@@H]1O)c1ccccc1 NAME: '
               '(1R,2S)-3-phenylcyclohexa-3,5-diene-1,2-diol REASON: MISSED '
               'Molecule contains 0 valid hydroxyl groups, diols must have at '
               'least two\n'
               ' * SMILES: OCCSCCO NAME: thiodiglycol REASON: MISSED Molecule '
               'contains 0 valid hydroxyl groups, diols must have at least '
               'two\n'
               ' * SMILES: OCC(O)COC(=O)c1ccccc1Nc1ccnc2cc(Cl)ccc12 NAME: '
               'glafenine REASON: MISSED Molecule contains 0 valid hydroxyl '
               'groups, diols must have at least two\n'
               ' * SMILES: OC[C@H](O)COc1ccc(Cl)cc1 NAME: (S)-chlorphenesin '
               'REASON: MISSED Molecule contains 0 valid hydroxyl groups, '
               'diols must have at least two\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)CO NAME: '
               '1,2-octacosanediol REASON: MISSED Molecule contains 0 valid '
               'hydroxyl groups, diols must have at least two\n'
               ' * SMILES: '
               '[H][C@@]12OC(=O)C(=C)[C@@]1([H])[C@H](OC(=O)C(C)CC)C(=O)[C@@H](C)C[C@H](O)C[C@](C)(O)[C@H]2OC(=O)C(C)C '
               'NAME: cardivin B REASON: MISSED Molecule contains 0 valid '
               'hydroxyl groups, diols must have at least two\n'
               ' * SMILES: '
               '[H][C@]12CC(=C)[C@@H](O)CC\\C(C)=C\\[C@@H](O)[C@]1([H])C(=C)C(=O)O2 '
               'NAME: 1-epi-tatridin B REASON: MISSED Molecule contains 0 '
               'valid hydroxyl groups, diols must have at least two\n'
               ' * SMILES: '
               '[H][C@@]12[C@@H](O)[C@H](O)C(C)(C)[C@]11O[C@@](C)(O[C@@H]1CC1=CC=CC(C)=C1)[C@@H]2C '
               'NAME: mollebenzylanol B REASON: MISSED Molecule contains 0 '
               'valid hydroxyl groups, diols must have at least two\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C3=CC=C(C=C3)C(=O)N(C)C)O[C@@H]1CN(C)C(=O)NC4=CC=CC=C4F)[C@H](C)CO',
                                     'name': '4-[(4R,5S)-5-[[[(2-fluoroanilino)-oxomethyl]-methylamino]methyl]-2-[(2R)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-8-yl]-N,N-dimethylbenzamide',
                                     'reason': 'Molecule contains 1 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': 'C[C@H]1[C@H](CN(C(=O)C2=C(C=CC(=C2)C#N)OC[C@H]3[C@H](CC[C@@H](O3)CCN(C1=O)C)OC)C)OC',
                                     'name': 'LSM-10564',
                                     'reason': 'Molecule contains 0 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': 'CCC(C(=O)NC1=NN=C(S1)SCC)(C(F)(F)F)C(F)(F)F',
                                     'name': 'N-[5-(ethylthio)-1,3,4-thiadiazol-2-yl]-2,2-bis(trifluoromethyl)butanamide',
                                     'reason': 'Molecule contains 0 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': 'O1C(C(CC1=O)C)CCCCCC',
                                     'name': '5-Hexyldihydro-4-methyl-2(3H)-furanone',
                                     'reason': 'Molecule contains 0 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': 'O=C1O[C@]23C(=O)N[C@H]([C@@H]2[C@H](C)C(=C[C@@H]3C=CC[C@@H]([C@H]([C@H](C=C1)C)O)C)C)CC4=CC=CC=C4',
                                     'name': '7-deoxy-cytochalasin Z7',
                                     'reason': 'Molecule contains 1 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': '[C@]123[C@H](O1)[C@](C(C=C2[C@]4([C@](CC3)(C([C@@H](C(C4)=O)O)(C)C)[H])C)=O)(C=C)C',
                                     'name': 'yucalexin P8',
                                     'reason': 'Molecule contains 1 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': 'O=C1C2=C(C(OC)=CC=C2OC3=C1C(OC)=CC(=C3)C(=O)OC)C(=O)OC',
                                     'name': '2,8-dimethoxy-1,6-dimethoxycarbonyl-xanthen-9-one',
                                     'reason': 'Molecule contains 0 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': 'O([C@@H](C[N+](C)(C)C)CC([O-])=O)C(=O)CCCCCCCCCCCCCCCCCCC(O)=O',
                                     'name': 'CAR(DC20:0)',
                                     'reason': 'Molecule contains 1 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': 'C[C@H]1CN(C(=O)C2=CC=CC=C2C3=CC=CC=C3CO[C@@H]1CN(C)S(=O)(=O)C4=CC=CC(=C4)F)[C@H](C)CO',
                                     'name': 'LSM-12790',
                                     'reason': 'Molecule contains 1 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'},
                                 {   'smiles': 'CCCC#CC1=CC=C(C=C1)[C@H]2[C@@H](N([C@@H]2C#N)S(=O)(=O)C)CO',
                                     'name': '(2S,3R,4R)-4-(hydroxymethyl)-1-methylsulfonyl-3-(4-pent-1-ynylphenyl)azetidine-2-carbonitrile',
                                     'reason': 'Molecule contains 1 hydroxyl '
                                               'group(s), diols must have at '
                                               'least two'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12CCCN1C(=O)[C@]1(O)[C@@H](OCC=C(C)C)c3c4n([C@H](OOC(C)(C)C[C@]4([H])N1C2=O)C=C(C)C)c1cc(OC)ccc31',
                                      'name': 'fumitremorgin A',
                                      'reason': 'Molecule contains 1 hydroxyl '
                                                'group(s), diols must have at '
                                                'least two'}],
    'attempt': 2,
    'reasoning': '**Analysis of the previous program and outcomes:**\n'
                 '\n'
                 'The previous code attempts to identify hydroxyl groups (–OH) '
                 'by iterating over oxygen atoms and checking if they are '
                 'bonded to both a hydrogen and a carbon atom, where the '
                 'carbon is not part of a carbonyl group (C=O). However, RDKit '
                 'typically handles hydrogens implicitly unless they are '
                 'explicitly added using `Chem.AddHs()`. This means that '
                 '`neighbor.GetAtomicNum() == 1` will not find any hydrogens '
                 'because they are not present as explicit atoms in the '
                 'molecule.\n'
                 '\n'
                 'As a result, the `has_hydrogen` variable remains `False`, '
                 'and the code fails to recognize any hydroxyl groups, leading '
                 'to a hydroxyl count of zero for all molecules. This explains '
                 'why all the test molecules were incorrectly classified as '
                 'not being diols.\n'
                 '\n'
                 '**How to improve:**\n'
                 '\n'
                 'To accurately identify hydroxyl groups without relying on '
                 'explicit hydrogens, we can use SMARTS patterns that account '
                 'for implicit hydrogens. The SMARTS pattern `[OX2H]` matches '
                 'oxygen atoms that have two connected atoms (implicit or '
                 'explicit) and one implicit hydrogen, which corresponds to '
                 'hydroxyl groups. By using `mol.GetSubstructMatches()` with '
                 'this pattern, we can find all hydroxyl groups in the '
                 'molecule.\n'
                 '\n'
                 'We can then count the number of hydroxyl groups by counting '
                 'the number of matches to the pattern. If the molecule '
                 'contains at least two hydroxyl groups, we classify it as a '
                 'diol.\n'
                 '\n'
                 '**Updated code:**',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 400,
    'num_false_positives': 100,
    'num_true_negatives': 166,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.8,
    'recall': 0.9975062344139651,
    'f1': 0.8879023307436182,
    'accuracy': 0.848575712143928,
    'negative_predictive_value': 0.9940119760479041}