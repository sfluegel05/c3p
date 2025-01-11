"""
Classifies: CHEBI:35819 branched-chain fatty acid
"""
from rdkit import Chem

def is_branched_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a branched-chain fatty acid based on its SMILES string.
    A branched-chain fatty acid has one or more alkyl substituents and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a branched-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (R-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for evidence of branching by checking for non-linear shapes
    # Consider a simple search for carbon with more than two bonds as potential branch points
    branch_pattern = Chem.MolFromSmarts("[C]([C])([C])")
    if not mol.HasSubstructMatch(branch_pattern) and not Chem.MolFromSmarts("[CX4](C)C"):
        return False, "No branch points found"

    # Check for ring structures, considering that some BCFAs can indeed have small rings like cyclopropyl
    # We will allow BCFAs to have one or two rings as this is not uncommon in complex natural BCFAs
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 2:  
        return False, "Too many ring structures for a typical branched-chain fatty acid"

    return True, "Contains a carbon chain with branch points and a carboxylic acid group with acceptable complexity"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35819',
                          'name': 'branched-chain fatty acid',
                          'definition': 'Any fatty acid in which the parent '
                                        'hydrocarbon chain has one or more '
                                        'alkyl substituents; a common '
                                        'component in animal and bacterial '
                                        'lipids. The fatty acyl chain is '
                                        'usually saturated and the substituent '
                                        'a methyl group; however, unsaturated '
                                        'BCFAs are found in marine animals, '
                                        'and branches other than methyl are '
                                        'found in microbial lipids.',
                          'parents': ['CHEBI:35366'],
                          'xrefs': ['KEGG:C05996', 'PMID:18318842'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.5689655172413793 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: CC(C)CCCC(C)CCCC(C)CCCC(C)C(O)=O NAME: '
               'pristanic acid REASON: CORRECT Contains a carbon chain with '
               'branch points and a carboxylic acid group with acceptable '
               'complexity\n'
               ' * SMILES: CC(C)(C)CC(O)=O NAME: 3,3-dimethylbutyric acid '
               'REASON: CORRECT Contains a carbon chain with branch points and '
               'a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: OC(C(CC)C)C(O)=O NAME: 2-hydroxy-3-methylpentanoic '
               'acid REASON: CORRECT Contains a carbon chain with branch '
               'points and a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: CC(CO)CCCCCCCCCCCCCC(O)=O NAME: '
               'omega-hydroxy-15-methylpalmitic acid REASON: CORRECT Contains '
               'a carbon chain with branch points and a carboxylic acid group '
               'with acceptable complexity\n'
               ' * SMILES: CC(CCCCCC(=O)O)CC(=O)O NAME: 3-Methylazelaic acid '
               'REASON: CORRECT Contains a carbon chain with branch points and '
               'a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: CC(C)CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O NAME: '
               '26-methylheptacosanoic acid REASON: CORRECT Contains a carbon '
               'chain with branch points and a carboxylic acid group with '
               'acceptable complexity\n'
               ' * SMILES: CC(C)CCCCCCCCC(O)=O NAME: 10-methylundecanoic acid '
               'REASON: CORRECT Contains a carbon chain with branch points and '
               'a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: C([C@@H](CCC[C@H](CCC[C@H](CCCC(C)C)C)C)C)(=O)O '
               'NAME: (2R,6S,10S)-2,6,10,14-pristanic acid REASON: CORRECT '
               'Contains a carbon chain with branch points and a carboxylic '
               'acid group with acceptable complexity\n'
               ' * SMILES: CC(CC(C)C(=O)O)CC(=O)O NAME: 2,4-Dimethyladipic '
               'acid REASON: CORRECT Contains a carbon chain with branch '
               'points and a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCC(O)=O NAME: '
               '26-methyloctacosanoic acid REASON: CORRECT Contains a carbon '
               'chain with branch points and a carboxylic acid group with '
               'acceptable complexity\n'
               ' * SMILES: CCCC[C@H](C)[C@@H](O)CC(O)=O NAME: '
               '(3S,4S)-3-hydroxy-4-methyloctanoic acid REASON: CORRECT '
               'Contains a carbon chain with branch points and a carboxylic '
               'acid group with acceptable complexity\n'
               ' * SMILES: OC(=O)C(CCCC(C)C)C NAME: 2,6-dimethylheptanoic acid '
               'REASON: CORRECT Contains a carbon chain with branch points and '
               'a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: CC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O NAME: '
               '28-methylnonacosanoic acid REASON: CORRECT Contains a carbon '
               'chain with branch points and a carboxylic acid group with '
               'acceptable complexity\n'
               ' * SMILES: OC(=O)\\C=C\\C(C)(C)C NAME: '
               '4,4-dimethyl-2E-pentenoic acid REASON: CORRECT Contains a '
               'carbon chain with branch points and a carboxylic acid group '
               'with acceptable complexity\n'
               ' * SMILES: OC(=O)/C=C\\C(C)(C)C NAME: '
               '4,4-dimethyl-2Z-pentenoic acid REASON: CORRECT Contains a '
               'carbon chain with branch points and a carboxylic acid group '
               'with acceptable complexity\n'
               ' * SMILES: CCC(C)CCCCCCCCCCCCCCCCC(O)=O NAME: '
               '18-methylicosanoic acid REASON: CORRECT Contains a carbon '
               'chain with branch points and a carboxylic acid group with '
               'acceptable complexity\n'
               ' * SMILES: OC(=O)CCC/C=C\\CC/C=C\\C/C=C\\CCC(CC)C NAME: '
               '16-methyl-octadeca-5Z,9Z,12Z-trienoic acid REASON: CORRECT '
               'Contains a carbon chain with branch points and a carboxylic '
               'acid group with acceptable complexity\n'
               ' * SMILES: CC(C)CCCCCCCCCCCC(O)=O NAME: isopentadecanoic acid '
               'REASON: CORRECT Contains a carbon chain with branch points and '
               'a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: CCCCC(C)CCCCCC(C)CCCCCCCCCC(C)CC(O)=O NAME: '
               '3,13,19-trimethyltricosanoic acid REASON: CORRECT Contains a '
               'carbon chain with branch points and a carboxylic acid group '
               'with acceptable complexity\n'
               ' * SMILES: CC(C)CCCCCCCCCCCCCCC(O)=O NAME: '
               '16-methylheptadecanoic acid REASON: CORRECT Contains a carbon '
               'chain with branch points and a carboxylic acid group with '
               'acceptable complexity\n'
               ' * SMILES: OC(=O)C(CCC)(C)C NAME: alpha,alpha-dimethyl valeric '
               'acid REASON: CORRECT Contains a carbon chain with branch '
               'points and a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: CC(C)C[C@@H](O)C(O)=O NAME: '
               '(R)-2-hydroxy-4-methylpentanoic acid REASON: CORRECT Contains '
               'a carbon chain with branch points and a carboxylic acid group '
               'with acceptable complexity\n'
               ' * SMILES: CC(CCCCC(=O)O)CC(=O)O NAME: 3-Methylsuberic acid '
               'REASON: CORRECT Contains a carbon chain with branch points and '
               'a carboxylic acid group with acceptable complexity\n'
               ' * SMILES: CC(C)CCCCCCCCCCCCCCCCCCCCCCC(O)=O NAME: '
               '24-methylpentacosanoic acid REASON: CORRECT Contains a carbon '
               'chain with branch points and a carboxylic acid group with '
               'acceptable complexity\n'
               ' * SMILES: CC(C)C(C)C(O)=O NAME: 2,3-dimethylbutyric acid '
               'REASON: CORRECT Contains a carbon chain with branch points and '
               'a carboxylic acid group with acceptable complexity\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCC[C@H]([C@H](O)CCCCCCCCCCCCCCCCC[C@H]1C[C@H]1CCCCCCCCCCCCCCCC[C@@H](OC)[C@H](C)CCCCCCCCCCCCCCCCCC)C(O)=O '
               'NAME: '
               '(2R)-2-[(1R)-1-hydroxy-18-{(1S,2R)-2-[(17R,18R)-17-methoxy-18-methylhexatriacontyl]cyclopropyl}octadecyl]hexacosanoic '
               'acid REASON: MISSED Structure too complex for a typical '
               'branched-chain fatty acid\n'
               ' * SMILES: OC(=O)C(CCC)CC NAME: alpha-ethyl valeric acid '
               'REASON: MISSED No branch points found\n'
               ' * SMILES: '
               'C\\C(CO)=C/CC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C(O)=O NAME: '
               '(2E,6E,10E,14E)-omega-hydroxygeranylgeranic acid REASON: '
               'MISSED No branch points found\n'
               ' * SMILES: CC(C)=CCC\\C(C)=C\\CC\\C(C)=C\\CC\\C(C)=C\\C(O)=O '
               'NAME: (2E,6E,10E)-geranylgeranic acid REASON: MISSED No branch '
               'points found\n'
               ' * SMILES: '
               'C(CCCCCCCCC1C(CCCCCCCCCCC2C(CCCCCCCCCCCCCCCCCCC[C@@H](O)[C@H](C(O)=O)CCCCCCCCCCCCCCCCCCCCCCCC)C2)C1)CCCCCCCCCC '
               'NAME: '
               '(2R)-2-[(1R)-1-hydroxy-20-{2-[10-(2-nonadecylcyclopropyl)decyl]cyclopropyl}icosyl]hexacosanoic '
               'acid REASON: MISSED No branch points found\n'
               ' * SMILES: CCCC(CCC)C(O)=O NAME: valproic acid REASON: MISSED '
               'No branch points found\n'
               ' * SMILES: '
               'C(CCCCCC1C(C1)CCCCCCCCCCCCC[C@@H](O)[C@H](C(O)=O)CCCCCCCCCCCCCCCCCCCCCCCC)CCCCCCCCC2C(C2)CCCCCCCCCCCCCCCCCCCC '
               'NAME: '
               '(2R)-2-[(1R)-1-hydroxy-14-{2-[14-(2-icosylcyclopropyl)tetradecyl]cyclopropyl}tetradecyl]hexacosanoic '
               'acid REASON: MISSED No branch points found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC\\C=C/CCCCCCCC\\C=C/CCCCCCCCCCCCCCCCCCC[C@@H](O)[C@@H](CCCCCCCC)C(O)=O '
               'NAME: '
               '(2R,3R,23Z,33Z)-3-hydroxy-2-octyldopentaconta-23,33-dienoic '
               'acid REASON: MISSED No branch points found\n'
               ' * SMILES: OC(=O)\\C=C\\C(C)C NAME: 4-Methyl-2-pentenoic acid '
               'REASON: MISSED No branch points found\n'
               ' * SMILES: OC(=O)CCCCC(=C)C(=C)C(O)=O NAME: 2,3-Methylene '
               'suberic acid REASON: MISSED No branch points found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'Too many ring structures for a '
                                               'typical branched-chain fatty '
                                               'acid'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'Too many ring structures for a '
                                               'typical branched-chain fatty '
                                               'acid'},
                                 {   'smiles': 'FC(F)(F)c1nc2c(Cl)c(Cl)c(Cl)c(Cl)c2[nH]1',
                                     'name': '4,5,6,7-Tetrachloro-2-trifluoromethylbenzimidazole',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'Cl.OC[C@H]1O[C@@H]2[C@@H](Oc3nc(=N)ccn23)[C@@H]1O',
                                     'name': 'ancitabine hydrochloride',
                                     'reason': 'No carboxylic acid group '
                                               'found'}],
    'sample_false_negatives': [],
    'attempt': 4,
    'reasoning': '### Analysis of Previous Code and Outcomes\n'
                 'The previous attempt to classify branched-chain fatty acids '
                 '(BCFAs) using certain substructure searches and checks '
                 'failed, specifically because it did not adequately capture '
                 'the complexity needed to identify all examples correctly. '
                 'Here are key observations and potential improvements:\n'
                 '\n'
                 '1. **Branch Detection:**\n'
                 '   - The SMARTS pattern for branch detection seems to have '
                 'overlooked some structures. Alkyl side chains or other '
                 'substituents might not strictly correspond to the pattern '
                 'specified.\n'
                 '   - BCFAs are characterized by at least one branch, '
                 'typically a methyl group, but the placements and types of '
                 'branches can vary.\n'
                 '\n'
                 '2. **Ring Structures:**\n'
                 '   - The check for the allowable number of rings might be '
                 'too restrictive. Some complex BCFAs are naturally occurring, '
                 'which can include one or more cyclopropyl groups or other '
                 'small-ring structures.\n'
                 '\n'
                 '3. **Carbon Counts and Complexity:**\n'
                 '   - The threshold for carbon count and branch complexity '
                 'might need reevaluation to ensure it captures BCFAs '
                 'accurately.\n'
                 '   - Large and complex fatty acids might present with high '
                 'numbers of carbon atoms and branches which need to be '
                 'reflected in the logic.\n'
                 '\n'
                 'By refining the substructure patterns and conditions, '
                 'focusing more on the key features of BCFAs—alkyl '
                 'substituents and carboxylic groups—we can improve the '
                 'classification.\n'
                 '\n'
                 '### Improved Code\n'
                 'Here is an updated approach to classify BCFAs more '
                 'accurately:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 169,
    'num_false_positives': 100,
    'num_true_negatives': 325,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.6282527881040892,
    'recall': 1.0,
    'f1': 0.771689497716895,
    'accuracy': 0.8316498316498316,
    'negative_predictive_value': 1.0}