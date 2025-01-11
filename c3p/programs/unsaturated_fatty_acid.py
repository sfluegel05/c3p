"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    An unsaturated fatty acid contains at least one C=C or C#C bond and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an unsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "Missing carboxylic acid group"

    # Look for carbon-carbon double bond (C=C) or triple bond (C#C)
    unsaturation_pattern_1 = Chem.MolFromSmarts("C=C")
    unsaturation_pattern_2 = Chem.MolFromSmarts("C#C")
    
    if not mol.HasSubstructMatch(unsaturation_pattern_1) and not mol.HasSubstructMatch(unsaturation_pattern_2):
        return False, "No C=C or C#C bond found"

    return True, "Contains carboxylic acid group and C=C or C#C bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:27208',
                          'name': 'unsaturated fatty acid',
                          'definition': 'Any fatty acid containing at least '
                                        'one C=C or C#C bond.',
                          'parents': ['CHEBI:35366'],
                          'xrefs': [   'LIPID_MAPS_class:LMFA0103',
                                       'PMID:5322381'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No C=C or C#C bond found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'Missing carboxylic acid group'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No C=C or C#C bond found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'Missing carboxylic acid group'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'Missing carboxylic acid group'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'No C=C or C#C bond found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No C=C or C#C bond found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No C=C or C#C bond found'},
                                 {   'smiles': 'O1[C@H]([C@H](NC(=O)C)[C@@H](O)C[C@]1(O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O)C(O)=O)[C@H](O)[C@H](O[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](OC(=O)C)[C@H](O)CO)C(O)=O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-6-[(1S,2R)-2-[(2S,4S,5R,6R)-5-acetamido-6-[(1R,2R)-1-acetyloxy-2,3-dihydroxypropyl]-2-carboxy-4-hydroxyoxan-2-yl]oxy-1,3-dihydroxypropyl]-2-[(2R,3S,4S,5R,6S)-3,5-dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-4-hydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No C=C or C#C bond found'},
                                 {   'smiles': 'OCCCCCCCCC1CC2C3C(C4C3CC4)C2CC1',
                                     'name': '8-[3]-ladderane-1-octanol',
                                     'reason': 'Missing carboxylic acid '
                                               'group'}],
    'sample_false_negatives': [   {   'smiles': 'O=C1[C@@H]([C@@H](CC1)CC(O)=O)CCCCC',
                                      'name': '(+)-7-epi--9,10-dihydrojasmonic '
                                              'acid',
                                      'reason': 'No C=C or C#C bond found'},
                                  {   'smiles': 'O=C1[C@H]([C@H](CC1)CC(O)=O)CCCCC',
                                      'name': '(-)-7-epi--9,10-dihydrojasmonic '
                                              'acid',
                                      'reason': 'No C=C or C#C bond found'},
                                  {   'smiles': 'S(C(=O)CCCCC[C@@H]1[C@@H](C(C=C1)=O)C/C=C\\CC)CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP([O-])([O-])=O)O)(=O)[O-])(=O)[O-])(C)C)O)=O',
                                      'name': '(9S,13S)-1a,1b-dinor-12-oxo-10,15-phytodienoyl-CoA(4-)',
                                      'reason': 'Missing carboxylic acid '
                                                'group'},
                                  {   'smiles': 'O=C1N(C(=O)C=2N=CN=NC2N1C)C',
                                      'name': 'Fervenulin',
                                      'reason': 'Missing carboxylic acid '
                                                'group'},
                                  {   'smiles': 'S(C(CCCCCCC[C@@H]1[C@@H](C(C=C1)=O)C/C=C\\CC)=O)CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP([O-])([O-])=O)O)(=O)[O-])(=O)[O-])(C)C)O)=O',
                                      'name': '(15Z)-12-oxophyto-10,15-dienoyl-CoA(4-)',
                                      'reason': 'Missing carboxylic acid '
                                                'group'},
                                  {   'smiles': 'C12=C(N=CN=C1N)N(C=N2)[C@@H]3O[C@@H]([C@@H](C3O)OP(=O)([O-])[O-])COP(=O)(OP(=O)(OCC([C@](C(NCCC(NCCSC(CCCCCCC[C@@H]4[C@@H](C(CC4)=O)C/C=C\\CC)=O)=O)=O)(O)[H])(C)C)[O-])[O-]',
                                      'name': 'OPC8-CoA (4-)',
                                      'reason': 'Missing carboxylic acid '
                                                'group'},
                                  {   'smiles': 'C(C(=O)OC)[C@@H]1[C@H](C(CC1)=O)CCCCC',
                                      'name': 'Methyl dihydrojasmonate',
                                      'reason': 'No C=C or C#C bond found'},
                                  {   'smiles': 'OC(CCC1(CC)C(=O)NC(=O)NC1=O)(C)C',
                                      'name': "3'-Hydroxyamobarbital",
                                      'reason': 'Missing carboxylic acid '
                                                'group'},
                                  {   'smiles': 'O=C1C(C(CC1)CC(O)=O)CCCCC',
                                      'name': '9,10-Dihydrojasmonic acid',
                                      'reason': 'No C=C or C#C bond found'},
                                  {   'smiles': 'S(C(=O)CCCCC[C@@H]1[C@@H](C(CC1)=O)C/C=C\\CC)CCNC(CCNC(=O)[C@@H](C(COP(OP(OC[C@H]2O[C@@H](N3C4=C(C(=NC=N4)N)N=C3)[C@@H]([C@@H]2OP([O-])([O-])=O)O)(=O)[O-])(=O)[O-])(C)C)O)=O',
                                      'name': '(9S,13S)-1a,1b-dinor-10,11-dihydro-12-oxo-15-phytoenoyl-CoA(4-)',
                                      'reason': 'Missing carboxylic acid '
                                                'group'}],
    'attempt': 0,
    'reasoning': 'To classify chemical entities as unsaturated fatty acids '
                 'based on their SMILES strings, we need to identify specific '
                 'patterns in the molecule that indicate the presence of '
                 'unsaturation and a fatty acid structure.\n'
                 '\n'
                 '### Step-by-step Approach:\n'
                 '1. **Parse the SMILES String**: Convert the SMILES string '
                 'into a molecular structure using RDKit.\n'
                 '2. **Identify Fatty Acid Structure**: Ensure the presence of '
                 'a carboxylic acid group `(C(=O)O)`. This is the hallmark of '
                 'a fatty acid.\n'
                 '3. **Check for Unsaturation**: Look for carbon-carbon double '
                 'bonds `(C=C)` or triple bonds `(C#C)`. This signifies that '
                 'the fatty acid is unsaturated.\n'
                 '4. **Return the Classification**: If both the fatty acid '
                 'structure and unsaturation are present, classify it as an '
                 'unsaturated fatty acid.\n'
                 '\n'
                 'Here is the Python code implementing this logic:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 833,
    'num_false_positives': 100,
    'num_true_negatives': 542,
    'num_false_negatives': 10,
    'num_negatives': None,
    'precision': 0.8928188638799571,
    'recall': 0.9881376037959668,
    'f1': 0.938063063063063,
    'accuracy': 0.9259259259259259,
    'negative_predictive_value': 0.9818840579710145}