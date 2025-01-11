"""
Classifies: CHEBI:62499 methyl-branched fatty acid
"""
from rdkit import Chem

def is_methyl_branched_fatty_acid(smiles: str):
    """
    Determines if a molecule is a methyl-branched fatty acid based on its SMILES string.
    A methyl-branched fatty acid must have a carboxylic acid group and only methyl branches.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a methyl-branched fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the carboxylic acid group (must be terminal, C(=O)O must not have further connections)
    terminal_carboxylic_pattern = Chem.MolFromSmarts("C(=O)[O;H1]")
    if not mol.HasSubstructMatch(terminal_carboxylic_pattern):
        return False, "No terminal carboxylic acid group found"
    
    # Look for non-methyl branching
    for atom in mol.GetAtoms():
        # If a carbon atom has more than one non-hydrogen neighbor, it could be a branching point
        if atom.GetSymbol() == 'C' and atom.GetDegree() > 3:
            # Count how many of the neighbors are methyl-type carbons (CH3)
            methyl_branches = sum(1 for neighbor in atom.GetNeighbors() if neighbor.GetSymbol() == 'C' and neighbor.GetDegree() == 1)
            # If there is a branch that is not a methyl group, return False
            if atom.GetDegree() - methyl_branches > 2:
                return False, "Non-methyl branches detected"

    return True, "Contains only methyl branches with a terminal carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:62499',
                          'name': 'methyl-branched fatty acid',
                          'definition': 'Any branched-chain fatty acid '
                                        'containing methyl branches only.',
                          'parents': ['CHEBI:35819'],
                          'xrefs': [   'PMID:17030019',
                                       'PMID:19747846',
                                       'PMID:19933331'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: SMILES: O(C(CC(CC(CC(CCC)C)C)C)C)C=O NAME: '
               '1,3,5,7-Tetramethyldecyl formate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCC)C=O NAME: Hexyl formate REASON: WRONGLY '
               'CLASSIFIED Contains methyl branches only with a carboxylic '
               'acid group\n'
               ' * SMILES: O(CCCCCC/C=C\\CC/C=C\\CCCC)C(=O)C NAME: '
               '7Z,11Z-Hexadecadienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CC/C=C\\CCCCCCCC)C(=O)C NAME: 3Z-Dodecenyl '
               'acetate REASON: WRONGLY CLASSIFIED Contains methyl branches '
               'only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCC/C=C/C=C\\C=C)C(=O)C NAME: '
               '7E,9Z,11-Dodecatrienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CCC/C=C/C=C\\CC/C=C\\CCCCC)C(=O)C NAME: '
               '4E,6Z,10Z-Hexadecatrienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCC#C/C=C\\CC)C(=O)C NAME: 9Z-Dodecen-7-ynyl '
               'acetate REASON: WRONGLY CLASSIFIED Contains methyl branches '
               'only with a carboxylic acid group\n'
               ' * SMILES: O(C\\C=C\\CCC)C=O NAME: 2-Hexenyl formate REASON: '
               'WRONGLY CLASSIFIED Contains methyl branches only with a '
               'carboxylic acid group\n'
               ' * SMILES: O(CCCCC/C=C/CC)C(=O)C NAME: 6E-Nonenyl acetate '
               'REASON: WRONGLY CLASSIFIED Contains methyl branches only with '
               'a carboxylic acid group\n'
               ' * SMILES: C(=C\\COC(C)=O)(\\CCC=C(C)C)/C NAME: neryl acetate '
               'REASON: WRONGLY CLASSIFIED Contains methyl branches only with '
               'a carboxylic acid group\n'
               ' * SMILES: O(CC/C=C\\CCCCCCCC/C=C\\CCCC)C(=O)C NAME: '
               '3Z,13Z-Octadecadienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(C/C=C\\CC/C=C\\CC)C(=O)C NAME: '
               '(2E,6Z)-2,6-Nonadien-1-Yl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCCCCCCCCCC)C(=O)C NAME: Pentadecyl acetate '
               'REASON: WRONGLY CLASSIFIED Contains methyl branches only with '
               'a carboxylic acid group\n'
               ' * SMILES: O(CCCCCC\\C=C/C=C/C=C)C(=O)C NAME: '
               '7Z,9E,11-Dodecatrienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CC\\C=C\\CCCCCCCC/C=C\\CCCC)C(=O)C NAME: '
               '3E,13Z-Octadecadienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCC\\C=C/C=C\\C)C(=O)C NAME: '
               '8Z,10Z-Dodecadienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCCCCC/C=C\\C/C=C/C)C(=O)C NAME: '
               '11Z,14E-Hexadecadienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CC\\C=C\\CCCCCCCC)C(=O)C NAME: 3E-Dodecenyl '
               'acetate REASON: WRONGLY CLASSIFIED Contains methyl branches '
               'only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCC/C=C\\CCCCCC)C(=O)C NAME: 8Z-Pentadecenyl '
               'acetate REASON: WRONGLY CLASSIFIED Contains methyl branches '
               'only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCCCCC/C=C\\CCCCC)C(=O)C NAME: '
               '11Z-Heptadecenyl acetate REASON: WRONGLY CLASSIFIED Contains '
               'methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCCCC/C=C/C=C\\C)C(=O)C NAME: '
               '10E,12Z-Tetradecadienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               ' * SMILES: CCCCCCOC(C)=O NAME: hexyl acetate REASON: WRONGLY '
               'CLASSIFIED Contains methyl branches only with a carboxylic '
               'acid group\n'
               ' * SMILES: O(CCC\\C=C\\CCCCCCC)C(=O)C NAME: 4E-Dodecenyl '
               'acetate REASON: WRONGLY CLASSIFIED Contains methyl branches '
               'only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCCCCCCC/C=C/CCCC)C(=O)C NAME: '
               '13E-Octadecenyl acetate REASON: WRONGLY CLASSIFIED Contains '
               'methyl branches only with a carboxylic acid group\n'
               ' * SMILES: O(CCCCCCCCCC\\C=C/C=C\\CC)C(=O)C NAME: '
               '11Z,13Z-Hexadecadienyl acetate REASON: WRONGLY CLASSIFIED '
               'Contains methyl branches only with a carboxylic acid group\n'
               'False negatives: SMILES: CC(C)CCCC(C)CCCC(C)CCCC(C)C(O)=O '
               'NAME: pristanic acid REASON: MISSED Non-methyl branches '
               'detected\n'
               ' * SMILES: CC(C)(C)CC(O)=O NAME: 3,3-dimethylbutyric acid '
               'REASON: MISSED Non-methyl branches detected\n'
               ' * SMILES: CC(CO)CCCCCCCCCCCCCC(O)=O NAME: '
               'omega-hydroxy-15-methylpalmitic acid REASON: MISSED Non-methyl '
               'branches detected\n'
               ' * SMILES: OC(=O)C(CCC)CC NAME: alpha-ethyl valeric acid '
               'REASON: MISSED Non-methyl branches detected\n'
               ' * SMILES: CC(C)CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O NAME: '
               '26-methylheptacosanoic acid REASON: MISSED Non-methyl branches '
               'detected\n'
               ' * SMILES: CC(C)CCCCCCCCC(O)=O NAME: 10-methylundecanoic acid '
               'REASON: MISSED Non-methyl branches detected\n'
               ' * SMILES: '
               'C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C(O)=O NAME: '
               '(2E,6E,10E,14E)-omega-hydroxygeranylgeranic acid REASON: '
               'MISSED Non-methyl branches detected\n'
               ' * SMILES: CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C(O)=O '
               'NAME: (2E,6E,10E)-geranylgeranic acid REASON: MISSED '
               'Non-methyl branches detected\n'
               ' * SMILES: C([C@@H](CCC[C@H](CCC[C@H](CCCC(C)C)C)C)C)(=O)O '
               'NAME: (2R,6S,10S)-2,6,10,14-pristanic acid REASON: MISSED '
               'Non-methyl branches detected\n'
               ' * SMILES: OC(=O)\\C=C\\C(C)C NAME: 4-Methyl-2-pentenoic acid '
               'REASON: MISSED Non-methyl branches detected\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': 'FC(F)(F)c1nc2c(Cl)c(Cl)c(Cl)c(Cl)c2[nH]1',
                                     'name': '4,5,6,7-Tetrachloro-2-trifluoromethylbenzimidazole',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'},
                                 {   'smiles': 'Cl.OC[C@H]1O[C@@H]2[C@@H](Oc3nc(=N)ccn23)[C@@H]1O',
                                     'name': 'ancitabine hydrochloride',
                                     'reason': 'No terminal carboxylic acid '
                                               'group found'}],
    'sample_false_negatives': [   {   'smiles': 'OC(=O)C(CC)(CC)C',
                                      'name': '2-ethyl-2-methyl-butanoic acid',
                                      'reason': 'Non-methyl branches detected'},
                                  {   'smiles': 'C(C(=CC(=O)[O-])C)CC=C(C)C',
                                      'name': '3,7-dimethylocta-2,6-dienoate(1-)',
                                      'reason': 'No terminal carboxylic acid '
                                                'group found'},
                                  {   'smiles': '[O-]C(=O)CCC([N+](C)(C)C)C',
                                      'name': '4-aminovaleric acid betaine',
                                      'reason': 'No terminal carboxylic acid '
                                                'group found'},
                                  {   'smiles': '[H][C@@]12[C@H](CCN1CC=C2COC(=O)[C@](O)([C@H](C)O)C(C)(C)O)OC(=O)C(\\C)=C/C',
                                      'name': 'heliosupine',
                                      'reason': 'No terminal carboxylic acid '
                                                'group found'}],
    'attempt': 1,
    'reasoning': 'To address the failures in classifying methyl-branched fatty '
                 'acids, we need to carefully revise the program with the goal '
                 'of differentiating these molecules accurately. The initial '
                 'approach correctly identified the presence of a carboxylic '
                 'acid group and checked for methyl branches, but it failed in '
                 'two aspects: identifying non-methyl branches wrongly and '
                 'classifying esters as methyl-branched fatty acids '
                 'mistakenly.\n'
                 '\n'
                 '### Hypotheses and Proposed Improvements:\n'
                 '1. **False Positives**: \n'
                 '   - The method classified esters as methyl-branched fatty '
                 'acids because having a carbonyl connected to an alkoxy group '
                 '(`O`) pattern was not properly distinguished from a '
                 'carboxylic acid group. This happened because of the broad '
                 'matching for carboxylic acid-like structures (e.g., '
                 'esters).\n'
                 '   - ***Improvement***: Ensure that the carbon in the '
                 'carboxylic group is terminal (or only has the carboxyl '
                 "branching), which doesn't typically leave the molecule open "
                 'for reclassification as esters.\n'
                 '\n'
                 '2. **True/False Negatives**: \n'
                 '   - Misclassifying structures without adequate '
                 'differentiation of backbone or placement of methyl branches '
                 'demonstrates there might be an oversight in detection. The '
                 'function seems overly restrictive when classifying '
                 'non-terminal carbons with methyl groups.\n'
                 '   - ***Improvement***: Check that the branches are indeed '
                 'exclusively methyl groups and ensure these groups are '
                 'accurately identified as branches from the main carbon '
                 'chain, not mistaken.\n'
                 '\n'
                 '3. **Implementation Details**:\n'
                 '   - Modify the program to look for the specific conditions '
                 'that define a "methyl-branched fatty acid," making sure the '
                 'branching carbons are attached directly to the backbone or '
                 'main chain, not forming secondary structures (like esters or '
                 'unrelated branching).\n'
                 '   - Improve the SMARTS patterns to correctly recognize '
                 'these features, and check correctly for the presence of '
                 'exactly one carboxylic acid end group (the terminal acid '
                 'functionality).\n'
                 '\n'
                 "Let's write the corrected program with these improvements:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 90,
    'num_false_positives': 100,
    'num_true_negatives': 553,
    'num_false_negatives': 4,
    'num_negatives': None,
    'precision': 0.47368421052631576,
    'recall': 0.9574468085106383,
    'f1': 0.6338028169014084,
    'accuracy': 0.8607764390896921,
    'negative_predictive_value': 0.992818671454219}