"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
from rdkit import Chem

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    A polyunsaturated fatty acid (PUFA) contains more than one carbon-carbon double bond and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyunsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group pattern (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"

    # Count the number of carbon-carbon double bonds (C=C)
    double_bond_count = len(mol.GetSubstructMatches(Chem.MolFromSmarts("C=C")))
    if double_bond_count <= 1:
        return False, f"Only {double_bond_count} double bond(s) found, need more than 1 for polyunsaturation"
    
    return True, "Contains more than one carbon-carbon double bond and a carboxylic acid group, indicative of a polyunsaturated fatty acid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26208',
                          'name': 'polyunsaturated fatty acid',
                          'definition': 'Any fatty acid containing more than '
                                        'one double bond. Acids in this group '
                                        'are reported to have cardioprotective '
                                        'effects; and levels are lowered in '
                                        'chronic fatigue syndrome.',
                          'parents': ['CHEBI:27208'],
                          'xrefs': [   'PMID:14977874',
                                       'PMID:16380690',
                                       'PMID:17891522'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'No carboxylic acid functional '
                                               'group found'},
                                 {   'smiles': 'CCCCCCCCCCCCCC(=O)OC[C@H](COP([O-])(=O)OC[C@@H](O)CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                     'name': '1-myristoyl-2-oleoyl-sn-glycero-3-phosphatidylglycerol(1-)',
                                     'reason': 'Only 1 double bond(s) found, '
                                               'need more than 1 for '
                                               'polyunsaturation'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'No carboxylic acid functional '
                                               'group found'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'No carboxylic acid functional '
                                               'group found'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'No carboxylic acid functional '
                                               'group found'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No carboxylic acid functional '
                                               'group found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No carboxylic acid functional '
                                               'group found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No carboxylic acid functional '
                                               'group found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'No carboxylic acid functional '
                                               'group found'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'No carboxylic acid functional '
                                               'group found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCC#CC\\C=C/CCCCCCCC(O)=O',
                                      'name': '(9Z)-octadec-9-en-12-ynoic acid',
                                      'reason': 'Only 1 double bond(s) found, '
                                                'need more than 1 for '
                                                'polyunsaturation'},
                                  {   'smiles': 'CCCCC\\C=C/CCC(O)=O',
                                      'name': 'cis-4-decenoic acid',
                                      'reason': 'Only 1 double bond(s) found, '
                                                'need more than 1 for '
                                                'polyunsaturation'},
                                  {   'smiles': 'OC(CCC1(CC)C(=O)NC(=O)NC1=O)(C)C',
                                      'name': "3'-Hydroxyamobarbital",
                                      'reason': 'No carboxylic acid functional '
                                                'group found'},
                                  {   'smiles': 'O=C1N(C(=O)C=2N=CN=NC2N1C)C',
                                      'name': 'Fervenulin',
                                      'reason': 'No carboxylic acid functional '
                                                'group found'},
                                  {   'smiles': 'OC(CCCCC)/C=C/C(=O)CCCCCCCCC(O)=O',
                                      'name': '(11E)-13-hydroxy-10-oxo-11-octadecenoic '
                                              'acid',
                                      'reason': 'Only 1 double bond(s) found, '
                                                'need more than 1 for '
                                                'polyunsaturation'},
                                  {   'smiles': 'CCCCCCCCC#CCC#CCC#CCCCC(O)=O',
                                      'name': '5,8,11-icosatriynoic acid',
                                      'reason': 'Only 0 double bond(s) found, '
                                                'need more than 1 for '
                                                'polyunsaturation'},
                                  {   'smiles': 'CCCC\\C=C/C#CC#CC(O)CCCCCCC(O)=O',
                                      'name': '(13Z)-8-hydroxyoctadecene-9,11-diynoic '
                                              'acid',
                                      'reason': 'Only 1 double bond(s) found, '
                                                'need more than 1 for '
                                                'polyunsaturation'},
                                  {   'smiles': 'C(CCCCCCCCC#CC#CCCCCC)(=O)O',
                                      'name': 'octadeca-10,12-diynoic acid',
                                      'reason': 'Only 0 double bond(s) found, '
                                                'need more than 1 for '
                                                'polyunsaturation'}],
    'attempt': 0,
    'reasoning': 'To classify a chemical entity as a polyunsaturated fatty '
                 'acid based on its SMILES string, the following '
                 'considerations must be taken into account:\n'
                 '\n'
                 '1. **Polyunsaturation**: The compound must have more than '
                 'one carbon-carbon double bond. This is the defining '
                 'characteristic of a polyunsaturated fatty acid.\n'
                 '\n'
                 '2. **Fatty Acid Structure**: It generally consists of a long '
                 'hydrocarbon chain with a terminal carboxylic acid group.\n'
                 '\n'
                 '3. **Examples**: The examples provided all have multiple '
                 'double bonds in their hydrocarbon chains along with a '
                 'carboxylic acid functional group.\n'
                 '\n'
                 'In the function, we will:\n'
                 '- Use RDKit to parse the SMILES string into a molecular '
                 'object.\n'
                 '- Check for the presence of a carboxylic acid group to '
                 'confirm the structure is a fatty acid.\n'
                 '- Count the number of double bonds in the molecule to '
                 'determine if it is polyunsaturated (more than one double '
                 'bond).\n'
                 '\n'
                 "Let's implement the function:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 600,
    'num_false_positives': 100,
    'num_true_negatives': 863,
    'num_false_negatives': 8,
    'num_negatives': None,
    'precision': 0.8571428571428571,
    'recall': 0.9868421052631579,
    'f1': 0.9174311926605504,
    'accuracy': 0.9312539783577339,
    'negative_predictive_value': 0.9908151549942594}