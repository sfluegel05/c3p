"""
Classifies: CHEBI:33709 amino acid
"""
from rdkit import Chem

def is_amino_acid(smiles: str):
    """
    Determines if a molecule is an amino acid based on its SMILES string.
    An amino acid is characterized by the presence of at least one carboxylic acid group 
    and one or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Look for amino group pattern
    amino_group_pattern = Chem.MolFromSmarts("N")
    if not mol.HasSubstructMatch(amino_group_pattern):
        return False, "No amino group found"
    
    return True, "Contains both carboxylic acid and amino groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33709',
                          'name': 'amino acid',
                          'definition': 'A carboxylic acid containing one or '
                                        'more amino groups.',
                          'parents': ['CHEBI:33575', 'CHEBI:50047'],
                          'xrefs': ['Wikipedia:Amino_acid'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@@H](O)CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                     'name': '1-myristoyl-2-oleoyl-sn-glycero-3-phosphatidylglycerol(1-)',
                                     'reason': 'No amino group found'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'No carboxylic acid group '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'CC(C)(S)[C@H](C(O)=O)n1ccnc1Cc1ccccc1',
                                      'name': 'benzylpenillamine',
                                      'reason': 'No amino group found'},
                                  {   'smiles': 'CC(=NNC(=O)C1=CC=C(C=C1)N(C)C)C2=CC(=CC=C2)OC',
                                      'name': '4-(dimethylamino)-N-[1-(3-methoxyphenyl)ethylideneamino]benzamide',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'C1=COC(=C1)NC(=O)C2=CC=C(C=C2)N',
                                      'name': '4-amino-N-(2-furanyl)benzamide',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'OC(=O)[C@@H]1CCn2c1ccc2C(=O)c1ccccc1',
                                      'name': '(R)-ketorolac',
                                      'reason': 'No amino group found'},
                                  {   'smiles': 'OC(=O)[C@H]1CCn2c1ccc2C(=O)c1ccccc1',
                                      'name': '(S)-ketorolac',
                                      'reason': 'No amino group found'},
                                  {   'smiles': 'CNC(=O)c1ccccc1N',
                                      'name': 'N-Methylanthranilamide',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'CC1=NC(=CC=C1)NC(=O)C2=CC(=C(C=C2)N(C)C3CCCCC3)[N+](=O)[O-]',
                                      'name': '4-[cyclohexyl(methyl)amino]-N-(6-methyl-2-pyridinyl)-3-nitrobenzamide',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'O=C1[C@H]2[C@@](CC[C@H]3[C@H]2[C@@H](O)[C@H](O)[C@H]4[C@@H](O)[C@H](O)CC[C@]34C)(C)[C@H](C1)[C@H]([C@@H]5OC(O)C[C@@H](C5)C(C)C)C',
                                      'name': 'Contignasterol',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'OC(=O)C1CCn2c1ccc2C(=O)c1ccccc1',
                                      'name': '5-benzoyl-2,3-dihydro-1H-pyrrolizine-1-carboxylic '
                                              'acid',
                                      'reason': 'No amino group found'},
                                  {   'smiles': 'C[n+]1ccccc1C([O-])=O',
                                      'name': 'Betaine homarine',
                                      'reason': 'No amino group found'}],
    'attempt': 0,
    'reasoning': 'To classify a chemical entity as an amino acid based on its '
                 'SMILES string, we need to identify two key structural '
                 'features: \n'
                 '\n'
                 '1. The presence of a carboxylic acid group, typically '
                 'represented as `C(=O)O` in SMILES notation.\n'
                 '2. The presence of at least one amino group, represented as '
                 '`N` in SMILES.\n'
                 '\n'
                 'Using RDKit, we can identify these substructures within the '
                 'given SMILES string. \n'
                 '\n'
                 'Here is the Python code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1094,
    'num_false_positives': 100,
    'num_true_negatives': 369,
    'num_false_negatives': 16,
    'num_negatives': None,
    'precision': 0.916247906197655,
    'recall': 0.9855855855855856,
    'f1': 0.9496527777777777,
    'accuracy': 0.9265357821405953,
    'negative_predictive_value': 0.9584415584415584}