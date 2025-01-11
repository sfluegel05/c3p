"""
Classifies: CHEBI:21575 N-acetyl-amino acid
"""
"""
Classifies: CHEBI:21547 N-acetyl-amino acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acetyl_amino_acid(smiles: str):
    """
    Determines if a molecule is an N-acetyl-amino acid based on its SMILES string.
    An N-acetyl-amino acid is an N-acyl-amino acid that has acetyl as the acyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an N-acetyl-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the acetyl group (CH3-CO-)
    acetyl_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])")
    if not mol.HasSubstructMatch(acetyl_pattern):
        return False, "No acetyl group found"

    # Look for a carboxyl group (COOH) which is part of the amino acid structure
    carboxyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H,OX1H0-]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No carboxyl group found"

    # Look for a nitrogen that is part of the amino acid structure
    nitrogen_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4]")
    if not mol.HasSubstructMatch(nitrogen_pattern):
        return False, "No nitrogen found in amino acid structure"

    # Check if the acetyl group is attached to a nitrogen that is part of the amino acid structure
    acetyl_nitrogen_pattern = Chem.MolFromSmarts("[CH3][CX3](=[OX1])[NX3;H2,H1,H0][CX4]")
    if not mol.HasSubstructMatch(acetyl_nitrogen_pattern):
        return False, "Acetyl group not attached to a nitrogen in the amino acid structure"

    return True, "Contains an acetyl group attached to a nitrogen in an amino acid structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:21575',
                          'name': 'N-acetyl-amino acid',
                          'definition': 'An N-acyl-amino acid that has acetyl '
                                        'as the acyl group.',
                          'parents': [   'CHEBI:22160',
                                         'CHEBI:22195',
                                         'CHEBI:51569'],
                          'xrefs': ['KEGG:C02847'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.44919786096256686 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: [C@H](CC(C)C)(NC(=O)C)C(O)=O NAME: '
               'N-acetyl-D-leucine REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: OC(=O)[C@@H](NC(=O)C)CC=1N(C=NC1)C NAME: '
               'N-acetyl-3-methyl-L-histidine REASON: CORRECT Contains an '
               'acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: CC(=O)N[13C@@H]([13CH2][13C](O)=O)[13C](O)=O NAME: '
               'N-acetyl-L-aspartic acid-1,2,3,4-(13)C4 REASON: CORRECT '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: C(CCNC(=O)C)[C@@H](C(O)=O)NC(=O)C NAME: bisorcic '
               'REASON: CORRECT Contains an acetyl group attached to the '
               'nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](Cc1cnc[nH]1)C(O)=O NAME: '
               'N-acetyl-L-histidine REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: [H]C(=O)CCC[C@H](NC(C)=O)C(O)=O NAME: '
               'N-acetyl-L-2-aminoadipic acid semialdehyde REASON: CORRECT '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: CC(=O)N1CCC[C@H]1C(O)=O NAME: N-acetyl-L-proline '
               'REASON: CORRECT Contains an acetyl group attached to the '
               'nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)NC(Cc1ccccc1)C(O)=O NAME: '
               'N-acetylphenylalanine REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)NC(Cc1c[nH]c2ccccc12)C(O)=O NAME: '
               'N-acetyltryptophan REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)NC(Cc1ccc(O)cc1)C(O)=O NAME: N-acetyltyrosine '
               'REASON: CORRECT Contains an acetyl group attached to the '
               'nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](CCCNC(N)=N)C(O)=O NAME: '
               'N(alpha)-acetyl-L-arginine REASON: CORRECT Contains an acetyl '
               'group attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: C[C@H](NC(C)=O)C(O)=O NAME: N-acetyl-L-alanine '
               'REASON: CORRECT Contains an acetyl group attached to the '
               'nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](Cc1ccccc1)C(O)=O NAME: '
               'N-acetyl-L-phenylalanine REASON: CORRECT Contains an acetyl '
               'group attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: O=C(O)[C@@H](NC(=O)C)C(C)C NAME: N-acetyl-L-valine '
               'REASON: CORRECT Contains an acetyl group attached to the '
               'nitrogen of an amino acid backbone\n'
               ' * SMILES: [C@@H](CC(C)C)(NC(=O)C)C(O)=O NAME: '
               'N-acetyl-L-leucine REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: C(O)([C@@H](NC(C)=O)CSC(C(O)=O)CC(=O)O)=O NAME: '
               'N-acetyl-S-(2-succino)-L-cysteine REASON: CORRECT Contains an '
               'acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: C(CC(C)C)(NC(=O)C)C(O)=O NAME: N-acetyl-leucine '
               'REASON: CORRECT Contains an acetyl group attached to the '
               'nitrogen of an amino acid backbone\n'
               ' * SMILES: CSCC[C@H](NC(C)=O)C(O)=O NAME: '
               'N-acetyl-L-methionine REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: '
               '[H]C(=O)c1ccc2oc3cc(=O)c(N)c(SC[C@H](NC(C)=O)C(O)=O)c3nc2c1 '
               'NAME: grixazone A REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](CCN)C(O)=O NAME: '
               '(2S)-2-acetamido-4-aminobutanoic acid REASON: CORRECT Contains '
               'an acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: CC(=O)N[C@H](Cc1ccccc1)C(O)=O NAME: '
               'N-acetyl-D-phenylalanine REASON: CORRECT Contains an acetyl '
               'group attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: C(=O)(C(CC=1NC=NC1)NC(=O)C)O NAME: '
               'N-acetylhistidine REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)NCC(O)=O NAME: N-acetylglycine REASON: '
               'CORRECT Contains an acetyl group attached to the nitrogen of '
               'an amino acid backbone\n'
               ' * SMILES: CSCC[C@@H](NC(C)=O)C(O)=O NAME: '
               'N-acetyl-D-methionine REASON: CORRECT Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](CC(O)=O)C(O)=O NAME: '
               'N-acetyl-L-aspartic acid REASON: CORRECT Contains an acetyl '
               'group attached to the nitrogen of an amino acid backbone\n'
               'False positives: SMILES: '
               'O([C@H]1[C@@H]([C@H]([C@H](O[C@H]2[C@@H]([C@H]([C@H](O)[C@H](O2)CO)O)NC(=O)C)[C@H](O1)CO)O)NC(=O)C)[C@@H]3[C@H](O[C@@H](O[C@@H]([C@@H]([C@H](C(=O)[O-])NC(=O)C)O)[C@H](O)CO)[C@@H]([C@H]3O)NC(=O)C)CO '
               "NAME: N,N',N'',N'''-tetraacetylchitotetraonate REASON: WRONGLY "
               'CLASSIFIED Contains an acetyl group attached to the nitrogen '
               'of an amino acid backbone\n'
               ' * SMILES: [S@@](C[C@@H](C([O-])=O)NC(=O)C)(CC=C)=O NAME: '
               'N-acetylalliin(1-) REASON: WRONGLY CLASSIFIED Contains an '
               'acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: COC1=CC(CC(NC(C)=O)C(O)=O)=CC=C1O NAME: '
               'N-Acetylvanilalanine REASON: WRONGLY CLASSIFIED Contains an '
               'acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: CC(=O)N[C@@H](CCCC[NH3+])C([O-])=O NAME: '
               'N(2)-acetyl-L-lysine zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: CNC(=O)SC[C@H](NC(C)=O)C(O)=O NAME: AMCC REASON: '
               'WRONGLY CLASSIFIED Contains an acetyl group attached to the '
               'nitrogen of an amino acid backbone\n'
               ' * SMILES: '
               'S(C(=O)[C@H]([C@@H](O)C(=O)OC)C(C)C)C[C@H](NC(=O)C)C(=O)O '
               'NAME: JBIR-69 REASON: WRONGLY CLASSIFIED Contains an acetyl '
               'group attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](CSCC(=O)c1ccccc1)C(O)=O NAME: '
               'N-acetyl-S-phenacyl-L-cysteine REASON: WRONGLY CLASSIFIED '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: OC(=O)[C@H](NC(=O)C)[C@H](CC)C NAME: '
               'N-acetyl-D-alloisoleucine REASON: WRONGLY CLASSIFIED Contains '
               'an acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: O=C([O-])[C@H](NC(=O)C)CCC(=O)[O-] NAME: '
               'N-acetyl-D-glutamate(2-) REASON: WRONGLY CLASSIFIED Contains '
               'an acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: S(C[C@H](NC(=O)C)C(O)=O)CCC#N NAME: '
               'N-Acetyl-S-(2-cyanoethyl)-cysteine REASON: WRONGLY CLASSIFIED '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: '
               '[C@H]([C@H]([C@@H]([C@H](/C=C/CCCCCC[C@@H](CCCCCC)O)O)O)O)(C(=O)[O-])NC(=O)C '
               'NAME: sphingofungin D zwitterion REASON: WRONGLY CLASSIFIED '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: S(C[C@H](NC(=O)C)C(O)=O)C(=S)NCC=C NAME: '
               'N-Acetyl-S-(N-allylthiocarbamoyl)cysteine REASON: WRONGLY '
               'CLASSIFIED Contains an acetyl group attached to the nitrogen '
               'of an amino acid backbone\n'
               ' * SMILES: '
               'S(C1=C2NC(=O)C=CC=CC=C[C@@H]([C@@H](O)CC(=O)C(C)=CC[C@@H](C=C[C@@H]([C@@H]([C@H](C=C(C(C=3C(C1=O)=C(C2=O)C=C(C)C3O)=O)C)C)O)C)O)C)CC(NC(=O)C)C(=O)O '
               'NAME: Diastovaricin II REASON: WRONGLY CLASSIFIED Contains an '
               'acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: OC(=O)C(NC(=O)C)CC(O)=O NAME: N-Acetyl-DL-aspartic '
               'acid REASON: WRONGLY CLASSIFIED Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: [O-]C(=O)[C@@H](NC(=O)C)CC=1N=CN(C1)C NAME: '
               'N-acetyl-1-methyl-L-histidinate REASON: WRONGLY CLASSIFIED '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: '
               'O=[N+]([O-])C1=C2C(NC=C2C[C@@H](NC(=O)C)C(=O)O)=CC=C1 NAME: '
               'N-Acetyl-4-nitrotryptophan REASON: WRONGLY CLASSIFIED Contains '
               'an acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               ' * SMILES: O=C(O)[C@@H](NC(=O)C)CC1=CC(O)=C(O)C=C1 NAME: '
               'N-Acetyl-3,4-dihydroxy-L-phenylalanine REASON: WRONGLY '
               'CLASSIFIED Contains an acetyl group attached to the nitrogen '
               'of an amino acid backbone\n'
               ' * SMILES: CC(=O)NC(CCP(C)([O-])=O)C([O-])=O NAME: '
               'N-acetylphosphinothricin(2-) REASON: WRONGLY CLASSIFIED '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](C(O)=O)C(C)(C)SN=O NAME: '
               'S-nitroso-N-acetyl-D-penicillamine REASON: WRONGLY CLASSIFIED '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: '
               'S(C1=C(NC)C(=O)C=2C=C3C[C@@H](CC(=O)O)O[C@H](C3=C(C2C1=O)O)CCC(C)C)CC(NC(=O)C)C(=O)O '
               'NAME: '
               '2-acetamido-3-[[(1S,3S)-3-(carboxymethyl)-10-hydroxy-7-(methylamino)-1-(3-methylbutyl)-6,9-dioxo-3,4-dihydro-1H-benzo[g]isochromen-8-yl]sulanyl]propanoic '
               'acid REASON: WRONGLY CLASSIFIED Contains an acetyl group '
               'attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](CS[C@H]1[C@H](O)C=Cc2ccccc12)C(O)=O '
               'NAME: '
               '(1R)-N-Acetyl-L-cysteinyl-(2R)-hydroxy-1,2-dihydronaphthalene '
               'REASON: WRONGLY CLASSIFIED Contains an acetyl group attached '
               'to the nitrogen of an amino acid backbone\n'
               ' * SMILES: SC[C@H](N(C(=O)NCCCC(=O)O)C(=O)C)C(=O)O NAME: '
               'Gaburedin E REASON: WRONGLY CLASSIFIED Contains an acetyl '
               'group attached to the nitrogen of an amino acid backbone\n'
               ' * SMILES: CC(=O)N[C@@H](CCCC(=O)OP(O)(O)=O)C(O)=O NAME: '
               'N-acetyl-L-2-aminoadipic acid 6-phosphate REASON: WRONGLY '
               'CLASSIFIED Contains an acetyl group attached to the nitrogen '
               'of an amino acid backbone\n'
               ' * SMILES: FC=1C(CC(NC(=O)C)C(O)=O)=CC=CC1 NAME: '
               'N-acetyl-o-fluoro-dl-phenylalanine REASON: WRONGLY CLASSIFIED '
               'Contains an acetyl group attached to the nitrogen of an amino '
               'acid backbone\n'
               ' * SMILES: '
               'S(CC1=CC(OC)=C2C=3C(O)=C4C(=O)CCC(C4=C(C3CC2=C1)O)OC)CC(NC(=O)C)C(=O)O '
               'NAME: Cysfluoretin REASON: WRONGLY CLASSIFIED Contains an '
               'acetyl group attached to the nitrogen of an amino acid '
               'backbone\n'
               'False negatives: SMILES: CC(=O)NCCC(O)=O NAME: '
               'N-acetyl-beta-alanine REASON: MISSED No amino acid backbone '
               'found\n'
               ' * SMILES: CC(=O)NCCC[C@H](N)C(O)=O NAME: '
               'N(5)-acetyl-L-ornithine REASON: MISSED Acetyl group not '
               'attached to the nitrogen of the amino acid backbone\n'
               ' * SMILES: CC(=O)NCCCC[C@H](N)C(O)=O NAME: '
               'N(6)-acetyl-L-lysine REASON: MISSED Acetyl group not attached '
               'to the nitrogen of the amino acid backbone\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC[C@H]1O[C@H](C[C@@H]1O)N1C=NC2=C1N=CNC[C@H]2O',
                                     'name': 'pentostatin',
                                     'reason': 'No acetyl group found'},
                                 {   'smiles': 'O1C(CCC=2C1=CC(OC3OC(C(O)C(O)C3O)C(O)=O)=C(C2OC)C=4C(=O)C=5C(OC4)=CC(O)=C(O)C5)(C)C',
                                     'name': '6-{[6-(6,7-dihydroxy-4-oxo-4H-chromen-3-yl)-5-methoxy-2,2-dimethyl-3,4-dihydro-2H-1-benzopyran-7-yl]oxy}-3,4,5-trihydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No acetyl group found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@@H](OC)C=CC=C(C[C@@H](C)[C@@H]([C@@H]([C@@H]([C@@H](C=C(C=C1OC)C)C)O)C)O)C)[C@H]([C@@H](O)[C@@H]([C@@]2(O[C@H](/C=C/C)[C@@H](C)[C@@H](C2)OC3OC(C(O)C(C3)O)C)O)C)C',
                                     'name': 'Concanamycin D',
                                     'reason': 'No acetyl group found'},
                                 {   'smiles': 'O1[C@@H](O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]3CO)O)[C@H](NC(=O)C)[C@@H](O[C@@H]4O[C@@H]([C@H](O)[C@H](O)[C@H]4O)CO)[C@H](O[C@@H]5O[C@H]([C@@H](O)[C@@H](O)[C@@H]5O)C)[C@H]1CO',
                                     'name': 'N-[(3R,4R,5S,6R)-5-[(2S,3R,4S,5S,6R)-4-[(2S,3R,4R,5S,6R)-3-Acetamido-6-(hydroxymethyl)-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-5-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,4-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No carboxyl group found'},
                                 {   'smiles': 'O=C(O)[C@]1([C@H]2[C@@](OC=3C=C4C5=C(C=CC=C5)NC4=CC3CC2)(CC[C@@H]1O)C)C',
                                     'name': 'Oxiamycin',
                                     'reason': 'No acetyl group found'},
                                 {   'smiles': 'C1CCC(C1)CC#CC2=CC=C(C=C2)[C@@H]3[C@@H]4CN(CC(=O)N4[C@@H]3CO)C(=O)C5=CC=C(C=C5)F',
                                     'name': '(6R,7R,8S)-7-[4-(3-cyclopentylprop-1-ynyl)phenyl]-4-[(4-fluorophenyl)-oxomethyl]-8-(hydroxymethyl)-1,4-diazabicyclo[4.2.0]octan-2-one',
                                     'reason': 'No acetyl group found'},
                                 {   'smiles': 'S(=O)(C(SSCCC)CC)CCC',
                                     'name': 'Propyl 1-(propylsulfinyl)propyl '
                                             'disulfide',
                                     'reason': 'No acetyl group found'},
                                 {   'smiles': 'ClC=1C(=O)[C@@H]([C@@](O)(C/C=C\\CCCCC)C1)C[C@H](OC(=O)C)[C@@H](OC(=O)C)CCCC(OC)=O',
                                     'name': 'punaglandin 6',
                                     'reason': 'No carboxyl group found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCC[C@H](O)C(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)[C@@H](O)CCCCCCCCCCCCCC',
                                     'name': 'Ins-1-P-Cer(t18:0/2-OH-26:0)',
                                     'reason': 'No acetyl group found'},
                                 {   'smiles': 'CCC(C)(C)C(=O)OC1CC(C=C2C1[C@H]([C@H](C=C2)C)CC[C@@H]3CC(CC(=O)O3)O)C',
                                     'name': '2,2-dimethylbutanoic acid '
                                             '[(7S,8S)-8-[2-[(2R)-4-hydroxy-6-oxo-2-oxanyl]ethyl]-3,7-dimethyl-1,2,3,7,8,8a-hexahydronaphthalen-1-yl] '
                                             'ester',
                                     'reason': 'No acetyl group found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'The previous code had a few issues that led to false '
                 'positives and false negatives:\n'
                 '\n'
                 '1. **False Positives**: The code incorrectly classified '
                 'molecules that had an acetyl group attached to a nitrogen, '
                 'but not necessarily the nitrogen of an amino acid backbone. '
                 'This led to molecules like '
                 "N,N',N'',N'''-tetraacetylchitotetraonate being incorrectly "
                 'classified as N-acetyl-amino acids.\n'
                 '\n'
                 '2. **False Negatives**: The code missed some valid '
                 'N-acetyl-amino acids because it was too strict in requiring '
                 'the acetyl group to be directly attached to the nitrogen of '
                 'the amino acid backbone. For example, N-acetyl-beta-alanine '
                 "was missed because it doesn't have a traditional amino acid "
                 'backbone.\n'
                 '\n'
                 'To improve the code, we need to:\n'
                 '- Relax the requirement for the amino acid backbone to '
                 'include more variations.\n'
                 '- Ensure that the acetyl group is attached to a nitrogen '
                 'that is part of the amino acid structure, but not '
                 'necessarily the main backbone nitrogen.\n'
                 '- Consider the presence of a carboxyl group (COOH) as part '
                 'of the amino acid structure.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 45,
    'num_false_positives': 100,
    'num_true_negatives': 3894,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3103448275862069,
    'recall': 1.0,
    'f1': 0.4736842105263158,
    'accuracy': 0.9752413963852439,
    'negative_predictive_value': 1.0}