"""
Classifies: CHEBI:18035 diglyceride
"""
from rdkit import Chem

def is_diglyceride(smiles: str):
    """
    Determines if a molecule is a diglyceride based on its SMILES string.
    A diglyceride has a glycerol backbone with two ester-linked fatty acid chains.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with at least two oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("C(=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2 or len(ester_matches) > 2:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 2"

    return True, "Contains glycerol backbone with two ester-linked fatty acid chains"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18035',
                          'name': 'diglyceride',
                          'definition': 'A glyceride that is glycerol in which '
                                        'any two of the hydroxy groups have '
                                        'been acylated. In the structure '
                                        'shown, two of the R groups (positions '
                                        'not specified) are acyl groups while '
                                        'the remaining R group can be either H '
                                        'or an alkyl group.',
                          'parents': ['CHEBI:47778', 'CHEBI:76578'],
                          'xrefs': ['KEGG:C00165', 'LIPID_MAPS_class:LMGL0201'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1C2=C(OC(=C1)C)C3=C(OC)C=C(OC)C=C3C(=C2O)C4=C5OC(=CC(C5=C(O)C=6C4=CC(OC)=CC6OC)=O)C',
                                     'name': 'Isonigerone',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C(=C\\C/C=C\\CCCC(NC)=O)\\C/C=C\\C/C=C\\CCCCC',
                                     'name': 'N-methyl arachidonoyl amine',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'OC1(O)[C@]23N(CC1)C(N(O)[C@H]([C@@]3(N=C(N2)N)[H])CO)=N',
                                     'name': 'Decarbamoylneosaxitoxin',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'S([C@H]1N(C(=O)/C(=C\\C2=CC=C(OCC=C(C)C)C=C2)/N(C1=O)C)C)C',
                                     'name': 'Fusaperazine F',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'Oc1ccccc1I',
                                     'name': '2-iodophenol',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C1N(CC(=O)N[C@H](C(=O)O[C@H]([C@@H](C(NC(C=C1)=C)=O)C)C(CCCCCCCCCCCCCC)C)C(O)C(=O)N)C',
                                     'name': 'Rakicidin H',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O(C1=C(OC)C=C(C2=C(O)C(OC)=C(C3=CC=CC=C3)C=C2OC)C=C1)C/C=C(/CO)\\C',
                                     'name': 'Prenylterphenyllin F',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@@H](O)[C@@H]1OC[C@H]3O[C@@H](OC[C@H]4O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]4O)[C@H](O)[C@@H](O)[C@@H]3O)CO[C@@H]5O[C@@H]([C@@H](O)[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)[C@H]5O)CO[C@@H]7O[C@@H]([C@@H](O)[C@H](O)[C@H]7O)CO',
                                     'name': '(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5S,6R)-6-[[(2R,3R,4S,5R,6R)-6-[[(2R,3R,4S,5R,6R)-3,5-Dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxymethyl]-3,5-dihydroxy-4-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-2-yl]oxymethyl]-3,4,5-trihydroxyoxan-2-yl]oxymethyl]oxane-2,3,4,5-tetrol',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO[C@@H]3O[C@H]([C@@H](O)[C@@H](O)[C@@H]3O)C)O)[C@@H]4O[C@@H]([C@@H](O)[C@H](O[C@H]5O[C@@H]([C@@H](O)[C@H](O)[C@@H]5O[C@@H]6O[C@@H]([C@@H](O[C@@H]7O[C@@H]([C@H](O)[C@H](O)[C@H]7O)CO)[C@H](O)[C@@H]6NC(=O)C)CO)CO)[C@@H]4O)CO[C@H]8O[C@@H]([C@@H](O)[C@H](O)[C@@H]8O[C@@H]9O[C@@H]([C@@H](O[C@@H]%10O[C@@H]([C@H](O)[C@H](O)[C@H]%10O)CO)[C@H](O)[C@H]9NC(=O)C)CO)CO',
                                     'name': 'Gal2GlcNAc2Man3GlcNAcFucGlcNAc',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'CN1CC2(CCN(CC2)S(=O)(=O)C3=CN(C=N3)C)C4=C([C@@H]1CO)NC5=C4C=CC(=C5)OC',
                                     'name': "[(1R)-7-methoxy-2-methyl-1'-[(1-methyl-4-imidazolyl)sulfonyl]-1-spiro[3,9-dihydro-1H-pyrido[3,4-b]indole-4,4'-piperidine]yl]methanol",
                                     'reason': 'No glycerol backbone found'}],
    'sample_false_negatives': [   {   'smiles': 'CCCCCCCCCCCCCCCCOC[C@H](CO)OC(=O)CCCCCCCCCCC',
                                      'name': '1-O-hexadecyl-2-O-lauroyl-sn-glycerol',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'C(CCCCCCC/C=C\\CCCCCCCC)(=O)O[C@@H](COP(=O)([O-])[O-])CO/C=C\\CCCCCCCCCCCCCCCC',
                                      'name': '1-[(1Z)-octadecenyl]-2-oleoyl-sn-glycero-3-phosphate(2-)',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OCC(O/C=C\\CCCCCCCCCCCCCCCC)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)([O-])=O',
                                      'name': 'PC(22:6(4Z,7Z,10Z,13Z,16Z,19Z)/P-18:0)',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCOC[C@H](CO)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC',
                                      'name': '1-palmityl-2-arachidonoyl-sn-glycerol',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OCC(O/C=C\\CCCCCC/C=C\\CCCCCCCC)COC(=O)CCCCCCCCCCCCCC)([O-])=O',
                                      'name': 'PC(15:0/P-18:1(9Z))',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCOC[C@H](COP([O-])([O-])=O)OC(C)=O',
                                      'name': '1-palmityl-2-acetyl-sn-glycero-3-phosphate(2-)',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCOC[C@H](CO)OC(=O)CCCCCCC\\C=C/CCCCCCCC',
                                      'name': '1-palmityl-2-oleoyl-sn-glycerol',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'P(OCC[N+](C)(C)C)(OCC(O/C=C\\CCCCCCCCCCCCCC)COC(=O)CCCCCCCCCCCCC)([O-])=O',
                                      'name': 'PC(14:0/P-16:0)',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCOC[C@H](COP([O-])([O-])=O)OC(=O)CCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC',
                                      'name': '1-palmityl-2-arachidonoyl-sn-glycero-3-phosphate(2-)',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'},
                                  {   'smiles': 'O(CCC([N+](C)(C)C)C([O-])=O)CC(O)COC(=O)CCCCCCC/C=C\\CCCCCCCC',
                                      'name': '4-[2-hydroxy-3-[(Z)-octadec-9-enoyl]oxypropoxy]-2-(trimethylazaniumyl)butanoate',
                                      'reason': 'Found 1 ester groups, need '
                                                'exactly 2'}],
    'attempt': 0,
    'reasoning': 'To classify a molecule as a diglyceride based on its SMILES '
                 'string, we need to follow certain key rules. A diglyceride '
                 'is defined as a glycerol backbone with two acyl groups '
                 '(fatty acid chains) attached, while the third hydroxyl group '
                 'might be free or occupied by hydrogen or an alkyl group.\n'
                 '\n'
                 "1. **Parse the SMILES String:** We'll use RDKit to parse the "
                 'SMILES string into a molecular structure.\n'
                 '\n'
                 '2. **Identify Glycerol Backbone:** Look for a glycerol '
                 'backbone pattern, typically represented as [C-C-C] with each '
                 'carbon atom having at least one oxygen attached.\n'
                 '\n'
                 '3. **Identify Acyl Groups:** Look for two separate ester '
                 'groups, "-OC(=O)-" indicating an acyl chain bound to the '
                 'glycerol backbone. These groups are required for the '
                 'classification as a diglyceride.\n'
                 '\n'
                 '4. **Check for Third Group:** Ensure that the third hydroxyl '
                 'group remains either unmodified or modified by hydrogen or '
                 'an alkyl group.\n'
                 '\n'
                 '5. **Molecular Validation:** Ensure no more than two ester '
                 'groups are observed and classify accordingly.\n'
                 '\n'
                 "Here's how this pattern recognition and validation could be "
                 'implemented in a Python program:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 643,
    'num_false_positives': 100,
    'num_true_negatives': 3231,
    'num_false_negatives': 15,
    'num_negatives': None,
    'precision': 0.8654104979811574,
    'recall': 0.9772036474164134,
    'f1': 0.9179157744468237,
    'accuracy': 0.971170719478566,
    'negative_predictive_value': 0.9953789279112755}