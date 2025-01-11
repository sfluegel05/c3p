"""
Classifies: CHEBI:26004 phenylpropanoid
"""
"""
Classifies: phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    A phenylpropanoid is any organic aromatic compound with a structure based on a phenylpropane skeleton.
    This includes compounds like flavonoids, coumarins, lignins, stilbenes, and others with a C6-C3 backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define generalized SMARTS patterns for phenylpropanoid structures
    patterns = {
        'phenylpropane_core': Chem.MolFromSmarts('c1ccccc1CCC'),  # Phenyl ring linked to 3-carbon chain
        'phenylpropene_core': Chem.MolFromSmarts('c1ccccc1C=CC'),  # Phenyl ring linked to propenyl chain
        'phenylpropyne_core': Chem.MolFromSmarts('c1ccccc1C#CC'),  # Phenyl ring linked to propynyl chain
        'cinnamic_acid': Chem.MolFromSmarts('c1ccccc1C=CC(=O)O'),  # Cinnamic acid scaffold
        'coumarin_core': Chem.MolFromSmarts('O=C1C=CC2=CC=CC=C2O1'),  # Coumarin core
        'flavonoid_core': Chem.MolFromSmarts('c1cc(c(cc1)-c1coc2c1ccc(=O)c(=O)c2)O'),  # Flavonoid skeleton
        'isoflavonoid_core': Chem.MolFromSmarts('c1cc(c(cc1)O)-c1coc2c1ccc(=O)c(=O)c2'),  # Isoflavonoid skeleton
        'stilbene_core': Chem.MolFromSmarts('c1ccccc1C=Cc2ccccc2'),  # Stilbene scaffold
        'lignan_core': Chem.MolFromSmarts('c1cc(c(cc1)O)C[C@H](C2=CC=C(C=C2)O)O'),  # Lignan structure
    }

    # Check for matches to any of the patterns
    for name, pattern in patterns.items():
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name.replace('_', ' ')} substructure"

    # Check for C6-C3 backbone (phenyl attached to 3-carbons)
    phenyl = Chem.MolFromSmarts('c1ccccc1')
    if mol.HasSubstructMatch(phenyl):
        # Find phenyl ring atoms
        phenyl_match = mol.GetSubstructMatch(phenyl)
        phenyl_atoms = set(phenyl_match)

        # Search for connected 3-carbon chain
        for atom_idx in phenyl_atoms:
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() not in phenyl_atoms:
                    # Check if neighbor is connected to a 3-carbon chain
                    chain = Chem.PathToSubmol(mol, [atom_idx, neighbor.GetIdx()])
                    if chain.GetNumHeavyAtoms() >= 3:
                        return True, "Contains phenylpropane skeleton"

    return False, "Does not contain phenylpropanoid substructure"

__metadata__ = {   
    'chemical_class': {   
        'name': 'phenylpropanoid',
        'definition': 'Any organic aromatic compound with a structure based on a phenylpropane skeleton. The class includes naturally occurring phenylpropanoid esters, flavonoids, anthocyanins, coumarins and many small phenolic molecules as well as their semi-synthetic and synthetic analogues. Phenylpropanoids are also precursors of lignin.',
        'parents': []},
    'config': {   
        'llm_model_name': 'lbl/claude-sonnet',
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
    'attempt': 1,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26004',
                          'name': 'phenylpropanoid',
                          'definition': 'Any organic aromatic compound with a '
                                        'structure based on a phenylpropane '
                                        'skeleton. The class includes '
                                        'naturally occurring phenylpropanoid '
                                        'esters, flavonoids, anthocyanins, '
                                        'coumarins and many small phenolic '
                                        'molecules as well as their '
                                        'semi-synthetic and synthetic '
                                        'analogues. Phenylpropanoids are also '
                                        'precursors of lignin.',
                          'parents': ['CHEBI:33659'],
                          'xrefs': ['Wikipedia:Phenylpropanoid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.7844827586206896 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O1C(CC(=O)C=2C1=CC(OC(=O)C3=CC=CC=C3)=CC2O)C4=CC=CC=C4 NAME: '
               'Pinocembrin 7-O-benzoate REASON: CORRECT Contains '
               'phenylpropane core substructure\n'
               ' * SMILES: O1C(C(O)CC2=C1C=C(OC)C=C2OC)C3=CC(O)=C(OC)C=C3 '
               "NAME: 3,3'-Dihydroxy-4',5,7-trimethoxyflavan REASON: CORRECT "
               'Contains phenylpropane core substructure\n'
               ' * SMILES: '
               'O1[C@@H]([C@H](O)[C@H](C=2C1=CC(O)=CC2O)CC(O)=O)C3=CC=C(O)C=C3 '
               "NAME: 3'-Deoxydryopteric acid REASON: CORRECT Contains "
               'phenylpropane core substructure\n'
               ' * SMILES: '
               'O(CC(C(CC1=CC=2OCOC2C=C1)CO)CC3=CC=4OCOC4C=C3)C(=O)C NAME: '
               'Hemiariensin REASON: CORRECT Contains phenylpropane core '
               'substructure\n'
               ' * SMILES: '
               'COc1c(ccc(O)c1CC=C(C)C)[C@@H]1COc2c(CC=C(C)C)c(O)ccc2C1=O '
               'NAME: lespeflorin D1 REASON: CORRECT Contains phenylpropane '
               'core substructure\n'
               ' * SMILES: COc1ccc(\\C=C/C)cc1OC NAME: cis-isomethyleugenol '
               'REASON: CORRECT Contains phenylpropene core substructure\n'
               ' * SMILES: '
               'C=1(C(O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)COC(/C=C/C3=CC=C(O)C(=C3)O)=O)=CC4=C([O+]1)C=C(O)C=C4O)C5=CC(=C(O)C(=C5)OC)OC '
               'NAME: malvidin 3-O-{6-O-[(E)-caffeoyl]-beta-D-glucoside} '
               'REASON: CORRECT Contains phenylpropene core substructure\n'
               ' * SMILES: COc1ccc(cc1)[C@@H]1CC(=O)c2c(O)cc(O)cc2O1 NAME: '
               "4'-methoxy-5,7-dihydroxyflavanone REASON: CORRECT Contains "
               'phenylpropane core substructure\n'
               ' * SMILES: '
               'O1C(C(O)C(C=2C=3OC(C(O)CC3C(O)=CC2O)C4=CC(O)=C(O)C=C4)C=5C1=CC(O)=CC5O)C6=CC(O)=C(O)C(O)=C6 '
               'NAME: Epigallocatechin-(4beta->8)-catechin REASON: CORRECT '
               'Contains phenylpropane core substructure\n'
               ' * SMILES: '
               'O=C1OC2=C(C=CC(=C2C)O[C@@H]3OC([C@H](OC)[C@H]([C@H]3O)OC(=O)N)(C)C)C(=C1NC(=O)C4=CC5=C(OC(C)(C)C(C5)O)C=C4)O '
               'NAME: Coumabiocin D REASON: CORRECT Contains phenylpropane '
               'core substructure\n'
               ' * SMILES: O1C(C(O)CC2=C1C=C(O)C=C2)C3=CC(O)=C(O)C=C3 NAME: '
               "3,3',4',7-Tetrahydroxyflavan REASON: CORRECT Contains "
               'phenylpropane core substructure\n'
               ' * SMILES: '
               'COc1ccc(cc1O)[C@@H]1CC(=O)c2c(O)cc(O[C@@H]3O[C@H](CO[C@@H]4O[C@@H](C)[C@H](O)[C@@H](O)[C@H]4O)[C@@H](O)[C@H](O)[C@H]3O)cc2O1 '
               'NAME: hesperidin REASON: CORRECT Contains phenylpropane core '
               'substructure\n'
               ' * SMILES: COc1cc(\\C=C/C(O)=O)cc(OC)c1O NAME: cis-sinapic '
               'acid REASON: CORRECT Contains phenylpropene core substructure\n'
               ' * SMILES: '
               'S(OCC1C(CC=2C(C1C3=CC(OC)=C(O)C(OC)=C3)=C(OC)C(O)=C(OC)C2)CO)(O)(=O)=O '
               "NAME: Lyoniresinol 9'-sulfate REASON: CORRECT Contains "
               'phenylpropane core substructure\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC2=CC=3OC=C(C4=C(OC)C=5CCC(OC5C=C4O)(C)C)C(=O)C3C=C2)C(O)=O '
               'NAME: '
               '3,4,5-trihydroxy-6-{[3-(7-hydroxy-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-6-yl)-4-oxo-4H-chromen-7-yl]oxy}oxane-2-carboxylic '
               'acid REASON: CORRECT Contains phenylpropane core substructure\n'
               ' * SMILES: '
               'CC(C)=CCC\\C(C)=C\\Cc1cc(ccc1O)[C@@H]1CC(=O)c2c(O1)cc(O)c(CC=C(C)C)c2O '
               'NAME: schizolaenone A REASON: CORRECT Contains phenylpropane '
               'core substructure\n'
               ' * SMILES: '
               'COc1ccc([C@@H]2COc3cc(O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4O)ccc3C2)c(O)c1OC '
               'NAME: '
               "(3R)-7,2'-dihydroxy-3',4'-dimethoxyisoflavan-7-O-beta-D-glucopyranoside "
               'REASON: CORRECT Contains phenylpropane core substructure\n'
               ' * SMILES: '
               'CCCCC[C@H]1CCCCCCCCCC(=O)O[C@H]2[C@@H](O)[C@@H](O[C@@H](C)[C@@H]2O[C@@H]2O[C@@H](C)[C@H](O[C@@H]3O[C@@H](C)[C@H](OC(=O)[C@@H](C)CC)[C@@H](OC(=O)\\C=C\\c4ccccc4)[C@H]3O)[C@@H](O[C@@H]3O[C@@H](C)[C@H](O)[C@@H](O)[C@H]3O)[C@H]2OC(=O)[C@@H](C)CC)O[C@@H]2[C@@H](O)[C@@H](O)[C@@H](C)O[C@H]2O1 '
               'NAME: pescaprein XXVII REASON: CORRECT Contains phenylpropene '
               'core substructure\n'
               ' * SMILES: OC1C(Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1 NAME: '
               'dihydromyricetin REASON: CORRECT Contains phenylpropane core '
               'substructure\n'
               ' * SMILES: '
               'O1C(CC(=O)C=2C1=CC(O)=CC2O)C3=CC(=C(OC)C(OC)=C3)CC=C(C)C NAME: '
               "5,7-Dihydroxy-3',4'-dimethoxy-5'-prenylflavanone REASON: "
               'CORRECT Contains phenylpropane core substructure\n'
               ' * SMILES: '
               'O1CC(CC=2C1=C(OC)C(OC)=C(OC)C2)C=3C(=O)C(OC)=C(OC)C(=O)C3 '
               "NAME: 6,7,8,3',4'-Pentamethoxyisoflavanquinone REASON: CORRECT "
               'Contains phenylpropane core substructure\n'
               ' * SMILES: '
               'O1C(=O)C([C@@H](CC(=O)C)C2=CC=C([N+]([O-])=O)C=C2)=C(O)C3=C1C=CC(O)=C3 '
               'NAME: 6-Hydroxy-R-acenocoumarol REASON: CORRECT Contains '
               'phenylpropane core substructure\n'
               ' * SMILES: COc1cc(cc(OC)c1O)\\C=C\\C(O)=O NAME: trans-sinapic '
               'acid REASON: CORRECT Contains phenylpropene core substructure\n'
               ' * SMILES: '
               'O[C@H]1O[C@H](COC(=O)\\C=C\\c2ccc(O)cc2)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: 6-O-(trans-4-coumaroyl)-alpha-D-glucopyranose REASON: '
               'CORRECT Contains phenylpropene core substructure\n'
               ' * SMILES: O(C(CCC=C(C)C)(C)C=C)C(=O)/C=C/C1=CC=CC=C1 NAME: '
               'Linalyl cinnamate REASON: CORRECT Contains phenylpropene core '
               'substructure\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'Oc1c(F)cc2cc(C(=O)ON3C(=O)CCC3=O)c(=O)oc2c1F NAME: pacific '
               'blue succinimidyl ester REASON: MISSED Does not contain '
               'phenylpropanoid substructure\n'
               ' * SMILES: O1C2=C(C(=O)C(C3=C(O)C(O)=C(OC)C=C3)=C1)C=CC(O)=C2 '
               'NAME: Koparin REASON: MISSED Does not contain phenylpropanoid '
               'substructure\n'
               ' * SMILES: CC(C)=CCc1c(O)ccc(c1O)-c1coc2cc(O)cc(O)c2c1=O NAME: '
               'licoisoflavone A REASON: MISSED Does not contain '
               'phenylpropanoid substructure\n'
               ' * SMILES: O1C=2C(=C(O)C=CC2)C=CC1=O NAME: '
               '5-hydroxy-2H-chromen-2-one REASON: MISSED Does not contain '
               'phenylpropanoid substructure\n'
               ' * SMILES: '
               'O1C(C(O)C(O)C(O)C1OC=2C=C3OC(=O)C(C(C)(C)C=C)=CC3=CC2OC)CO '
               'NAME: 3-(1,1-Dimethylallyl)scopoletin 7-glucoside REASON: '
               'MISSED Does not contain phenylpropanoid substructure\n'
               ' * SMILES: OC(C(O)CC/C(=C/COC1=CC=2OC(=O)C=CC2C=C1)/C)(C)C '
               'NAME: (R)-Marmin REASON: MISSED Does not contain '
               'phenylpropanoid substructure\n'
               ' * SMILES: '
               'O=C1OC2=C(C(=C(C3=C(O)C=C4OC(=O)C=C(C4=C3C)OC)C(=C2)O)C)C(=C1)OC '
               'NAME: Isokotanin C REASON: MISSED Does not contain '
               'phenylpropanoid substructure\n'
               ' * SMILES: COc1cc2oc(=O)c(cc2c(OC)c1CC=C(C)C)-c1ccc(O)cc1O '
               'NAME: glycyrin REASON: MISSED Does not contain phenylpropanoid '
               'substructure\n'
               ' * SMILES: '
               'C1=CC(=CC=C1C2=COC3=CC(=CC(=C3C2=O)O)OC4C(C(C(C(O4)CO)O)O)O)O '
               'NAME: '
               '5-hydroxy-3-(4-hydroxyphenyl)-7-[[3,4,5-trihydroxy-6-(hydroxymethyl)-2-oxanyl]oxy]-1-benzopyran-4-one '
               'REASON: MISSED Does not contain phenylpropanoid substructure\n'
               ' * SMILES: '
               'S(OC1=CC=C(C=2C(=O)C=3C(OC2)=CC(O)=CC3O)C=C1)(O)(=O)=O NAME: '
               '[4-(5,7-dihydroxy-4-oxo-4H-chromen-3-yl)phenyl]oxidanesulfonic '
               'acid REASON: MISSED Does not contain phenylpropanoid '
               'substructure\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C(CC)(O)(C)C',
                                     'name': '2-methylbutan-2-ol',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': 'O1[C@@H]([C@H](O)[C@@H](O)[C@@H](NC(=O)C)C1O)CO',
                                     'name': 'N-Acetyl-D-Gulosamine',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': '[C@@]12([C@]3([C@](CC[C@@]1([C@@]4(C(C[C@H](CC4)OS([O-])(=O)=O)=CC2)C)[H])([C@](CC3)([C@H](C)CC[C@@H](C(C)C)O)[H])C)[H])[H]',
                                     'name': '(24S)-hydroxycholesterol '
                                             '3-sulfate(1-)',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(C)O',
                                     'name': '17,20-dihydroxypregn-4-en-3-one',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': 'P(=O)(N=C(N(CC)CC)C)(OCC)F',
                                     'name': 'A-234 nerve agent',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': 'C[C@H](CC\\C=C(/C)CO[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@@]2(C)[C@@H]3CC=C4[C@@H](CC[C@H](O)C4(C)C)[C@]3(C)C(=O)C[C@]12C',
                                     'name': 'Carnosifloside I',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': 'OS(=O)(=O)CCS',
                                     'name': 'coenzyme M',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\C/C=C\\CCCCC)COC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\CCCCC)(O)(O)=O',
                                     'name': 'PA(20:3(8Z,11Z,14Z)/20:2(11Z,14Z))',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': 'O=C1O[C@H](CCC(=O)C[C@H](CC(=O)/C=C/C)C(O[C@H]([C@@H](C=C1)O)C)=O)C',
                                     'name': 'Acaulone A',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)COC(=O)CCCCCCCCCCC)(O)(O)=O',
                                     'name': 'PA(12:0/22:1(11Z))',
                                     'reason': 'Does not contain '
                                               'phenylpropanoid substructure'}],
    'sample_false_negatives': [   {   'smiles': 'O1[C@@]2([C@H](O)C[C@H]3C(CCC[C@@]3([C@H]2C[C@@]4([C@H]1[C@H](O)C(=C)[C@@H]([C@H]4O)O)O)C)(C)C)C',
                                      'name': 'Coniothyrenol',
                                      'reason': 'Does not contain '
                                                'phenylpropanoid '
                                                'substructure'}],
    'attempt': 2,
    'reasoning': "It appears that the previous code's SMARTS patterns were too "
                 'specific and did not capture the full diversity of the '
                 'phenylpropanoid class. Phenylpropanoids encompass a wide '
                 'variety of structures based on a C6-C3 phenylpropane '
                 'skeleton, including flavonoids, coumarins, isoflavones, '
                 'stilbenes, and more. The overly specific patterns failed to '
                 'match some valid phenylpropanoids, leading to false '
                 'negatives and a lower F1 score.\n'
                 '\n'
                 'To improve the classification, we can:\n'
                 '\n'
                 '- **Use more general SMARTS patterns** that capture the core '
                 'phenylpropanoid skeleton in its various forms, including '
                 'substitutions and ring systems.\n'
                 '- **Include additional patterns** for different subclasses '
                 'like isoflavonoids, coumarins, and stilbenes.\n'
                 '- **Allow flexibility** in the connectivity of the C6-C3 '
                 'unit to account for cyclizations and rearrangements common '
                 'in phenylpropanoids.\n'
                 '\n'
                 'By broadening the patterns and increasing their generality, '
                 'we can better identify phenylpropanoids and reduce false '
                 'negatives.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 2440,
    'num_false_positives': 100,
    'num_true_negatives': 89,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.9606299212598425,
    'recall': 0.9995903318312167,
    'f1': 0.9797229471993575,
    'accuracy': 0.961596958174905,
    'negative_predictive_value': 0.9888888888888889}