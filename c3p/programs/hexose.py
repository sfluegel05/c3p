"""
Classifies: CHEBI:18133 hexose
"""
"""
Classifies: CHEBI:24895 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.AllChem import ReactionFromSmarts

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde group (aldohexose)
    or a ketone group (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define hexose core patterns
    # Pyranose patterns (both alpha and beta)
    pyranose_patterns = [
        "[C]-1-[C]-[C]-[C]-[C](-[O]-1)(-[OH])-[CH2][OH]",  # Basic pyranose
        "[C]-1-[C]-[C]-[C]-[C](-[O]-1)-[CH2][OH]",         # Alternative form
        "[C]-1-[C]-[C]-[C]-[C](-[O]-1)-[CH2]-[O]"          # Modified form
    ]
    
    # Furanose patterns
    furanose_patterns = [
        "[C]-1-[C]-[C]-[C](-[O]-1)-[C](-[OH])-[CH2][OH]",  # Basic furanose
        "[C]-1-[C]-[C]-[C](-[O]-1)-[C]-[CH2][OH]",         # Alternative form
        "[C]-1-[C]-[C]-[C](-[O]-1)-[C]-[CH2]-[O]"          # Modified form
    ]
    
    # Open chain patterns (aldehexoses and ketohexoses)
    open_chain_patterns = [
        "[CH](=O)-[C]-[C]-[C]-[C]-[CH2][OH]",              # Aldehexose
        "[C](=O)-[C]-[C]-[C]-[C]-[CH2][OH]",               # Ketohexose
        "[CH](=O)-[C]-[C]-[C]-[C]-[CH2]-[O]"               # Modified form
    ]

    # Check for hexose core structure
    found_core = False
    for pattern in pyranose_patterns + furanose_patterns + open_chain_patterns:
        patt = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(patt):
            found_core = True
            break
    
    if not found_core:
        return False, "No hexose core structure found"

    # Count carbons in the core structure (should be 6)
    core_carbons = len(mol.GetSubstructMatch(patt))
    if core_carbons < 6:
        return False, f"Core structure contains fewer than 6 carbons ({core_carbons})"

    # Check for characteristic hydroxyl pattern
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    if mol.HasSubstructMatch(hydroxyl_pattern):
        hydroxyl_count = len(mol.GetSubstructMatches(hydroxyl_pattern))
        if hydroxyl_count < 2:  # Allow for modified hexoses
            return False, f"Too few hydroxyl groups ({hydroxyl_count})"
    
    # Check it's not a disaccharide or larger
    glycosidic_pattern = Chem.MolFromSmarts("[OX2]([CH1][OX2])[CH1]")
    if mol.HasSubstructMatch(glycosidic_pattern):
        glycosidic_matches = mol.GetSubstructMatches(glycosidic_pattern)
        if len(glycosidic_matches) > 1:
            return False, "Appears to be a polysaccharide"

    # Additional check for characteristic carbon chain
    chain_pattern = Chem.MolFromSmarts("[C]-[C]-[C]-[C]-[C]-[C]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic six-carbon chain"

    return True, "Contains hexose core structure with appropriate substitution pattern"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18133',
                          'name': 'hexose',
                          'definition': 'Any six-carbon monosaccharide which '
                                        'in its linear form contains either an '
                                        'aldehyde group at position 1 '
                                        '(aldohexose) or a ketone group at '
                                        'position 2 (ketohexose).',
                          'parents': ['CHEBI:35381'],
                          'xrefs': [   'GlyGen:G70994MS',
                                       'GlyTouCan:G70994MS',
                                       'KEGG:C00738'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.6437246963562753 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: OC[C@H]1OC(O)(CO)[C@H](O)[C@H]1O NAME: '
               'D-sorbofuranose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@]1(O)OC[C@@H](O)[C@H](O)[C@H]1O NAME: '
               'alpha-D-sorbopyranose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: O1[C@@H]([C@@H](O)[C@@H](O)[C@@H]1O)[C@@H](O)CO '
               'NAME: alpha-L-talofuranose REASON: CORRECT Matches hexose '
               'pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: [H]C(=O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)CO NAME: '
               'aldehydo-D-glucose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1O[C@@H](O)[C@H](O)C(=O)[C@H]1O NAME: '
               '3-keto-beta-D-galactose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: C1([C@H]([C@@H]([C@H]([C@H](O1)CO)O)O)O)O NAME: '
               'D-idopyranose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: O1[C@H]([C@@H](O)[C@H](O)[C@@H]1O)[C@@H](O)CO NAME: '
               'alpha-L-glucofuranose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H](O)[C@@H](O)[C@H](O)C(=O)CO NAME: '
               'keto-L-sorbose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H]1O NAME: '
               'alpha-D-tagatofuranose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: [C@H]1(O)C(O)O[C@H](C=O)[C@H]([C@@H]1O)O NAME: '
               '6-dehydro-D-glucose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O NAME: '
               'beta-D-allose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1O[C@](O)(CO)[C@H](O)[C@H]1O NAME: '
               'beta-D-sorbofuranose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: O1[C@@H]([C@@H](O)[C@H](O)[C@@H]1O)[C@H](O)CO NAME: '
               'beta-D-altrofuranose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1O[C@@H](O)[C@H](O)[C@H](O)[C@H]1O NAME: '
               'beta-D-gulose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@@H]1O[C@@H](O)[C@@H](O)[C@@H](O)[C@@H]1O NAME: '
               'alpha-L-gulose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]([C@H](CC(CO)=O)O)O NAME: '
               '3-deoxy-keto-D-fructose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: O1[C@@H]([C@H](O)[C@@H](O)[C@@H]1O)[C@@H](O)CO '
               'NAME: alpha-L-idofuranose REASON: CORRECT Matches hexose '
               'pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: [C@H](O)(C(CO)=O)[C@@H]([C@H](CO)O)O NAME: '
               'keto-L-tagatose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: O1[C@H]([C@@H](O)[C@H](O)C1O)[C@@H](O)CO NAME: '
               'L-glucofuranose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H](O)[C@@H](O)[C@@H](O)[C@@H](O)C=O NAME: '
               'aldehydo-L-talose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@@H]1OC(O)[C@@H](O)[C@H](O)[C@@H]1O NAME: '
               'L-galactopyranose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1O[C@](O)(CO)[C@@H](O)[C@H]1O NAME: '
               'beta-D-tagatofuranose REASON: CORRECT Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: [H]C(=O)[C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO NAME: '
               'aldehydo-L-gulose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: O1[C@H]([C@H](O)[C@H](O)[C@@H](O)C1O)CO NAME: '
               'L-altropyranose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               ' * SMILES: [H]C(=O)[C@H](O)[C@H](O)[C@@H](O)[C@H](O)CO NAME: '
               'aldehydo-D-gulose REASON: CORRECT Matches hexose pattern with '
               'correct number of carbons and hydroxyls\n'
               'False positives: SMILES: '
               '[NH3+][C@@H]1[C@@H](O)[C@@H](O)O[C@H](CO)[C@H]1O NAME: '
               'alpha-D-kanosamine(1+) REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: [H]C(=O)[C@H]([18F])[C@@H](O)[C@H](O)[C@H](O)CO '
               'NAME: 2-deoxy-2-((18)F)fluoro-aldehydo-D-glucose REASON: '
               'WRONGLY CLASSIFIED Matches hexose pattern with correct number '
               'of carbons and hydroxyls\n'
               ' * SMILES: O1[C@@H]([C@H](O)[C@H](O)[C@@H](O)C1O)C=O NAME: '
               'Galacto-hexodialdose REASON: WRONGLY CLASSIFIED Matches hexose '
               'pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@]1(O)[C@H](O)OC[C@@H](O)[C@H]1O NAME: '
               'beta-D-hamamelose REASON: WRONGLY CLASSIFIED Matches hexose '
               'pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: O1C(C(O)C(O)C(O)C1)CO NAME: 1,5-anhydrohexitol '
               'REASON: WRONGLY CLASSIFIED Matches hexose pattern with correct '
               'number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1OC(O)[C@H](F)[C@@H](O)[C@@H]1O NAME: '
               '2-deoxy-2-fluoro-D-glucopyranose REASON: WRONGLY CLASSIFIED '
               'Matches hexose pattern with correct number of carbons and '
               'hydroxyls\n'
               ' * SMILES: '
               'O([C@]([C@](O[H])([C@](O[H])(C(=O)[H])[H])[H])([C@@](O[H])(C([H])([H])[H])[H])[H])[H] '
               'NAME: 6-deoxy-L-allose REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: '
               'O([C@]([C@@](O[H])([C@](N([H])[H])(C(=O)[H])[H])[H])([C@@](O[H])(C(O[H])([H])[H])[H])[H])[H] '
               'NAME: L-glucosamine REASON: WRONGLY CLASSIFIED Matches hexose '
               'pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: '
               'O([C@@]([C@](O[H])(C(O[H])([H])[H])[H])([C@](O[H])(C(=O)C([H])([H])[H])[H])[H])[H] '
               'NAME: 1-deoxy-D-tagatose REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@@H](O)[C@H](O)CC(=O)C(O)=O NAME: '
               '2-dehydro-3-deoxy-D-galactonic acid REASON: WRONGLY CLASSIFIED '
               'Matches hexose pattern with correct number of carbons and '
               'hydroxyls\n'
               ' * SMILES: C(=O)(CC(=O)[C@H](O)[C@H](O)CO)O NAME: '
               '3-dehydro-2-deoxy-D-gluconic acid REASON: WRONGLY CLASSIFIED '
               'Matches hexose pattern with correct number of carbons and '
               'hydroxyls\n'
               ' * SMILES: O[C@@H](C[C@H](O)C(O)=O)[C@H](O)CO NAME: '
               '3-deoxy-D-arabino-hexonic acid REASON: WRONGLY CLASSIFIED '
               'Matches hexose pattern with correct number of carbons and '
               'hydroxyls\n'
               ' * SMILES: N[C@@H]1[C@@H](O)C(O)O[C@H](CO)[C@H]1O NAME: '
               '3-amino-3-deoxy-D-glucopyranose REASON: WRONGLY CLASSIFIED '
               'Matches hexose pattern with correct number of carbons and '
               'hydroxyls\n'
               ' * SMILES: [H]C(=O)[C@@](O)(CO)[C@H](O)[C@H](O)CO NAME: '
               'D-hamamelose REASON: WRONGLY CLASSIFIED Matches hexose pattern '
               'with correct number of carbons and hydroxyls\n'
               ' * SMILES: OCC(CO)(CO)NCC(O)=O NAME: '
               'N-tris(hydroxymethyl)methylglycine REASON: WRONGLY CLASSIFIED '
               'Matches hexose pattern with correct number of carbons and '
               'hydroxyls\n'
               ' * SMILES: [C@H]1(O)[C@H](O)[C@@H](CO)O[C@@H](O)[C@@H]1[NH3+] '
               'NAME: beta-D-glucosamine(1+) REASON: WRONGLY CLASSIFIED '
               'Matches hexose pattern with correct number of carbons and '
               'hydroxyls\n'
               ' * SMILES: N[C@H]([C@@H](O)[C@H](O)[C@H](O)CO)C([O-])=O NAME: '
               '2-amino-2-deoxy-D-gluconate REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: '
               'O1[C@]([C@@](O[H])([C@@](O[H])([C@](N([H])[H])(C1(O[H])[H])[H])[H])[H])(C(O[H])([H])[H])[H] '
               'NAME: L-galactosamine REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: N[C@@H](C=O)[C@@H](O)[C@@H](O)[C@H](O)CO NAME: '
               'aldehydo-D-galactosamine REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1OC[C@@H](O)[C@@H](O)[C@@H]1O NAME: '
               '1,5-anhydro-D-mannitol REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: NC(=O)[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@@H]1O '
               'NAME: beta-D-glucuronamide REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: C(=O)(C[C@@H](O)[C@H](O)[C@H](O)CO)[O-] NAME: '
               '2-deoxy-D-gluconate REASON: WRONGLY CLASSIFIED Matches hexose '
               'pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: '
               'O([C@@]([C@](O[H])([C@](O[H])(C(=O)[H])[H])[H])([C@@](O[H])(C([H])([H])[H])[H])[H])[H] '
               'NAME: 6-deoxy-L-gulose REASON: WRONGLY CLASSIFIED Matches '
               'hexose pattern with correct number of carbons and hydroxyls\n'
               ' * SMILES: OC[C@H]1O[C@H](O)[C@@H](F)[C@@H](O)[C@@H]1O NAME: '
               '2-deoxy-2-fluoro-alpha-D-mannose REASON: WRONGLY CLASSIFIED '
               'Matches hexose pattern with correct number of carbons and '
               'hydroxyls\n'
               ' * SMILES: [H][C@]1([18F])C(O)O[C@H](CO)[C@H](O)[C@@H]1O NAME: '
               '2-deoxy-2-((18)F)fluoro-D-galactopyranose REASON: WRONGLY '
               'CLASSIFIED Matches hexose pattern with correct number of '
               'carbons and hydroxyls\n'
               'False negatives: SMILES: '
               'O1C(OC)C(O)CC(C1C)(O)[C@@H](O)C[C@H](OC)[C@H](O)[C@H](O)C '
               'NAME: Caryophyllan REASON: MISSED Must have exactly 6 carbons, '
               'found 14\n'
               ' * SMILES: O[C@@H]([C@H](O)[C@H](O)C)[C@@H](O)CO NAME: '
               '6-Deoxyglucitol REASON: MISSED No sugar ring or carbonyl group '
               'found\n'
               ' * SMILES: '
               'O1C(O)[C@H](O)C[C@]([C@H]1C)(O)[C@@H](O)C[C@H](O)[C@H](O)[C@H](O)C '
               'NAME: Caryophyllose REASON: MISSED Must have exactly 6 '
               'carbons, found 12\n'
               ' * SMILES: '
               '[H]C([H])([C@@]([H])(O)C=O)[C@@]([H])(O)[C@@]([H])(C)O NAME: '
               'abequose REASON: MISSED Expected 5-6 oxygens, found 4\n'
               ' * SMILES: [H][C@](C)(O)[C@]([H])(O)[C@@]([H])(O)CC=O NAME: '
               'boivinose REASON: MISSED Expected 5-6 oxygens, found 4\n'
               ' * SMILES: O1C(C(O)C(O)C(O)C1OC(=O)/C(=C/CC=C(C)C)/C)CO NAME: '
               'Glucosyl (E)-2,6-Dimethyl-2,5-heptadienoate REASON: MISSED '
               'Must have exactly 6 carbons, found 15\n'
               ' * SMILES: CNCC(C(C(C(CO)O)O)O)O NAME: '
               '6-(methylamino)hexane-1,2,3,4,5-pentol REASON: MISSED Must '
               'have exactly 6 carbons, found 7\n'
               ' * SMILES: O1C(C(O)C(O)C(OC(=O)C2=CC=CC=C2)C1O)CO NAME: '
               '2-O-Benzoyl-D-glucose REASON: MISSED Must have exactly 6 '
               'carbons, found 13\n'
               ' * SMILES: C[C@@H](O)[C@@H](O)CCC=O NAME: amicetose REASON: '
               'MISSED Expected 5-6 oxygens, found 3\n'
               ' * SMILES: '
               'O=C(NCCC1=CC=CC=C1)[C@H](N[C@H](C(=O)NC=CC2=CC=C(O[C@@H]3O[C@@H]([C@H](NC(=O)C)[C@@H]([C@H]3O)O)C)C=C2)CC=4C5=C(C=CC=C5)N(C4)C)C(C)C '
               'NAME: Rhabdoplanin B REASON: MISSED Must have exactly 6 '
               'carbons, found 41\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'OC[C@H]1OC(O)(CO)[C@H](O)[C@H]1O',
                                      'name': 'D-sorbofuranose',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': 'O1C(OC)C(O)CC(C1C)(O)[C@@H](O)C[C@H](OC)[C@H](O)[C@H](O)C',
                                      'name': 'Caryophyllan',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': 'O[C@@H]([C@H](O)[C@H](O)C)[C@@H](O)CO',
                                      'name': '6-Deoxyglucitol',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': 'O1C(O)[C@H](O)C[C@]([C@H]1C)(O)[C@@H](O)C[C@H](O)[C@H](O)[C@H](O)C',
                                      'name': 'Caryophyllose',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': '[H]C([H])([C@@]([H])(O)C=O)[C@@]([H])(O)[C@@]([H])(C)O',
                                      'name': 'abequose',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@H](O)[C@@H](O)[C@H](O)C(=O)CO',
                                      'name': 'keto-L-sorbose',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@H]1O[C@@](O)(CO)[C@@H](O)[C@H]1O',
                                      'name': 'alpha-D-tagatofuranose',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': '[C@H]1(O)C(O)O[C@H](C=O)[C@H]([C@@H]1O)O',
                                      'name': '6-dehydro-D-glucose',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': '[H][C@](C)(O)[C@]([H])(O)[C@@]([H])(O)CC=O',
                                      'name': 'boivinose',
                                      'reason': 'No hexose core structure '
                                                'found'},
                                  {   'smiles': 'OC[C@H]1O[C@](O)(CO)[C@H](O)[C@H]1O',
                                      'name': 'beta-D-sorbofuranose',
                                      'reason': 'No hexose core structure '
                                                'found'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the issues with the previous code and suggest '
                 'improvements:\n'
                 '\n'
                 '1. **Carbon Count Issue**: The strict requirement of exactly '
                 '6 carbons is causing us to miss valid hexose derivatives and '
                 'modified hexoses. For example, "Caryophyllan" and '
                 '"2-O-Benzoyl-D-glucose" are valid hexoses with '
                 'modifications. We should focus on the core hexose structure '
                 'rather than total carbon count.\n'
                 '\n'
                 '2. **Oxygen Count**: The oxygen count check is too '
                 'restrictive. Modified hexoses may have different numbers of '
                 'oxygens due to substitutions or modifications.\n'
                 '\n'
                 '3. **Ring Pattern Detection**: The current ring pattern '
                 'detection is too simple. We need more specific patterns that '
                 'can identify the characteristic hexose ring structure with '
                 'its specific substitution pattern.\n'
                 '\n'
                 '4. **Derivatives Handling**: The code rejects molecules with '
                 'ester groups, but some valid hexose derivatives contain '
                 'esters or other modifications.\n'
                 '\n'
                 "Here's my improved approach:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 115,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.5348837209302325,
    'f1': 0.6969696969696969,
    'accuracy': 0.5348837209302325,
    'negative_predictive_value': 0.0}