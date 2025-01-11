"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is characterized as an aliphatic monocarboxylic acid with a chain of 4 to 28 carbons 
    (usually unbranched and even-numbered), which may be saturated or unsaturated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Get number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Fatty acids typically have between 4 and 28 carbon atoms
    # Adjust criteria to allow for more complexity than strictly chain length
    if c_count < 4:
        return False, f"Carbon chain length {c_count} too short for typical fatty acid"
    
    # We need to be cautious with long chains that might just be complex esters
    # Look at length but allow higher values with more checks
    if c_count > 28:
        larger_chain_info = f"Carbon chain length {c_count} for complexities beyond typical"
        # If structure is simple and exceeds typical complexity of these types
        if mol.GetRingInfo().NumRings() == 0:
            return False, larger_chain_info
        connected_heavy_atoms = mol.GetNumHeavyAtoms()
        if connected_heavy_atoms < c_count * 1.5: # Allow some connectivity above simple alkane chains
            return True, f"Larger but valid fatty acid: {larger_chain_info}"

    # Check for presence of rings and ensure they are small if present
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() > 0 and all(len(r) > 6 for r in ring_info.AtomRings()):
        return False, "Contains large ring structure(s), not characteristic of typical fatty acids"

    return True, "Valid fatty acid: Aliphatic monocarboxylic acid with primarily aliphatic character"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35366',
                          'name': 'fatty acid',
                          'definition': 'Any aliphatic monocarboxylic acid '
                                        'derived from or contained in '
                                        'esterified form in an animal or '
                                        'vegetable fat, oil or wax. Natural '
                                        'fatty acids commonly have a chain of '
                                        '4 to 28 carbons (usually unbranched '
                                        'and even-numbered), which may be '
                                        'saturated or unsaturated. By '
                                        'extension, the term is sometimes used '
                                        'to embrace all acyclic aliphatic '
                                        'carboxylic acids.',
                          'parents': ['CHEBI:18059', 'CHEBI:25384'],
                          'xrefs': [   'KEGG:C00162',
                                       'PMID:14287444',
                                       'PMID:14300208',
                                       'PMID:14328676',
                                       'Wikipedia:Fatty_acid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.8626373626373627 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: [H]C(CCCCCCCC)=C([H])CCCCCCCC(O)=O '
               'NAME: octadec-9-enoic acid REASON: CORRECT Valid fatty acid: '
               'Aliphatic monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: OC(=O)CCCCC\\C=C\\C=C(\\[N+]([O-])=O)/CCCCCC NAME: '
               'Dinor-10-NO2-CLA REASON: CORRECT Valid fatty acid: Aliphatic '
               'monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: OCC(C)=CC(O)=O NAME: 4-hydroxy-3-methylbut-2-enoic '
               'acid REASON: CORRECT Valid fatty acid: Aliphatic '
               'monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: CCCCCCCC[C@H](O)\\C=C\\CCCCCCC(O)=O NAME: '
               '(8E,10S)-10-hydroxy-8-octadecenoic acid REASON: CORRECT Valid '
               'fatty acid: Aliphatic monocarboxylic acid with between 4 and '
               '28 carbons\n'
               ' * SMILES: CCCCCCCCCCCC[C@@H](O)C(O)=O NAME: '
               '(2R)-2-hydroxytetradecanoic acid REASON: CORRECT Valid fatty '
               'acid: Aliphatic monocarboxylic acid with between 4 and 28 '
               'carbons\n'
               ' * SMILES: O[C@H](CCCCCCCC)/C=C/CCCCCCC(O)=O NAME: '
               '10R-HOME(8E) REASON: CORRECT Valid fatty acid: Aliphatic '
               'monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: OC(=O)CCCCCCCC=CCC=CCCCC NAME: '
               '9,12-Heptadecadienoic acid REASON: CORRECT Valid fatty acid: '
               'Aliphatic monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: OC(=O)CC[C@@H](CCCCCCCCCCCC)C NAME: '
               '4R-methyl-hexadecanoic acid REASON: CORRECT Valid fatty acid: '
               'Aliphatic monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: O=C(CCCC(O)=O)/C=C/C=C\\CCCCCCCCC NAME: '
               '(6E,8Z)-5-oxooctadecadienoic acid REASON: CORRECT Valid fatty '
               'acid: Aliphatic monocarboxylic acid with between 4 and 28 '
               'carbons\n'
               ' * SMILES: CCCCCCCCCCCCCCC(C)C(O)=O NAME: 2-methylhexadecanoic '
               'acid REASON: CORRECT Valid fatty acid: Aliphatic '
               'monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: OC(=O)C(CCCCCCCCCCCCCC(O)=O)C NAME: '
               '2-methyl-hexadecanedioic acid REASON: CORRECT Valid fatty '
               'acid: Aliphatic monocarboxylic acid with between 4 and 28 '
               'carbons\n'
               ' * SMILES: '
               'FC(F)(C(F)(F)C(F)(F)\\C(\\F)=C/C(O)=O)C(F)(F)C(F)(F)F NAME: '
               '(2E)-3,4,4,5,5,6,6,7,7,8,8,8-Dodecafluoro-2-octenoic acid '
               'REASON: CORRECT Valid fatty acid: Aliphatic monocarboxylic '
               'acid with between 4 and 28 carbons\n'
               ' * SMILES: OCCCCCCC(O)=O NAME: 7-hydroxyheptanoic acid REASON: '
               'CORRECT Valid fatty acid: Aliphatic monocarboxylic acid with '
               'between 4 and 28 carbons\n'
               ' * SMILES: OCCCCCCCCCCCCCCCCCCCCCC(O)=O NAME: '
               '22-hydroxydocosanoic acid REASON: CORRECT Valid fatty acid: '
               'Aliphatic monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: CCCC(CCC)C(O)=O NAME: valproic acid REASON: CORRECT '
               'Valid fatty acid: Aliphatic monocarboxylic acid with between 4 '
               'and 28 carbons\n'
               ' * SMILES: C[C@@H](O)CCCCCCCCCCCCCCCCCCC(O)=O NAME: '
               '(20R)-20-hydroxyhenicosanoic acid REASON: CORRECT Valid fatty '
               'acid: Aliphatic monocarboxylic acid with between 4 and 28 '
               'carbons\n'
               ' * SMILES: NCCCCCCCCCCC(O)=O NAME: 11-Aminoundecanoic acid '
               'REASON: CORRECT Valid fatty acid: Aliphatic monocarboxylic '
               'acid with between 4 and 28 carbons\n'
               ' * SMILES: CC\\C=C/C\\C=C/C\\C=C/CCCCCC[C@@H](O)C(O)=O NAME: '
               '(R)-2-hydroxy-alpha-linolenic acid REASON: CORRECT Valid fatty '
               'acid: Aliphatic monocarboxylic acid with between 4 and 28 '
               'carbons\n'
               ' * SMILES: CC[C@@H](C)C(O)=O NAME: (R)-2-methylbutyric acid '
               'REASON: CORRECT Valid fatty acid: Aliphatic monocarboxylic '
               'acid with between 4 and 28 carbons\n'
               ' * SMILES: [O-]C(=O)CCC/C=C\\CC NAME: cis-3-hexenylacetate '
               'REASON: CORRECT Valid fatty acid: Aliphatic monocarboxylic '
               'acid with between 4 and 28 carbons\n'
               ' * SMILES: OC(=O)CCCCCC#C NAME: 7-octynoic acid REASON: '
               'CORRECT Valid fatty acid: Aliphatic monocarboxylic acid with '
               'between 4 and 28 carbons\n'
               ' * SMILES: OC(=O)CCCCC/C=C/C\\C=C\\CCCCC NAME: '
               '7,10-hexadecadienoic acid REASON: CORRECT Valid fatty acid: '
               'Aliphatic monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: CC(=O)CCCC(O)=O NAME: 5-oxohexanoic acid REASON: '
               'CORRECT Valid fatty acid: Aliphatic monocarboxylic acid with '
               'between 4 and 28 carbons\n'
               ' * SMILES: OC(=O)CCCCCCCCCCCCC/C=C/CCCCCC NAME: 15-docosenoic '
               'acid REASON: CORRECT Valid fatty acid: Aliphatic '
               'monocarboxylic acid with between 4 and 28 carbons\n'
               ' * SMILES: OC(=O)CC\\C=C\\C\\C=C\\C\\C=C\\C\\C=C\\CCCCCCCC '
               'NAME: 4,7,10,13-Docosatetraenoic acid REASON: CORRECT Valid '
               'fatty acid: Aliphatic monocarboxylic acid with between 4 and '
               '28 carbons\n'
               'False positives: NONE\n'
               'False negatives: SMILES: O=C1OC(=C(C)C(=C1)OC)/C=C/C=C/C(=O)O '
               'NAME: Herbarin A REASON: MISSED Contains ring structure(s), '
               'usually not characteristic of typical fatty acids\n'
               ' * SMILES: '
               'C(CCC(O)=O)[C@@H](/C=C/C=C\\C/C=C\\C[C@@H]1[C@H](CCCCC)O1)OO '
               'NAME: '
               '(5S)-hydroperoxy-(14R,15S)-epoxy-(6E,8Z,11Z)-icosatrienoic '
               'acid REASON: MISSED Contains ring structure(s), usually not '
               'characteristic of typical fatty acids\n'
               ' * SMILES: O(C(CCCCCCCCCCCC(O)=O)CCC)C(=O)CCCCCCCCCCCCCCC '
               'NAME: FAHFA(16:0/13-O-16:0) REASON: MISSED Carbon chain length '
               '32 not in [4, 28]\n'
               ' * SMILES: O1C(CCCC)=CC=C1CCC2=CC(OC)=C(O)C=C2 NAME: '
               '2-Butyl-5-[2-(4-hydroxy-3-methoxyphenyl)ethyl]furan REASON: '
               'MISSED No carboxylic acid group found\n'
               ' * SMILES: '
               'CC\\C=C/CC1C(CCCCCCCC(=O)OC[C@H](CO[C@@H]2O[C@H](CO[C@H]3O[C@H](CO)[C@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H](O)[C@H]2O)OC(=O)CCCCCC2C=CC(=O)C2C\\C=C/CC)C=CC1=O '
               'NAME: arabidopside C REASON: MISSED Carbon chain length 49 not '
               'in [4, 28]\n'
               ' * SMILES: O=C1OC(CCCCC(CC(=CC(=CC(=O)O)C)C)C)C1CO NAME: 1233A '
               'REASON: MISSED Contains ring structure(s), usually not '
               'characteristic of typical fatty acids\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/[C@H](O)[C@H]1O[C@H]1CCCC(O)=O NAME: '
               '(7R)-hydroxy-(5S,6S)-epoxy-(8Z,11Z,14Z)-icosatrienoic acid '
               'REASON: MISSED Contains ring structure(s), usually not '
               'characteristic of typical fatty acids\n'
               ' * SMILES: CCCCC(C(=O)O)SC1=NC(=O)C(=C(N1)C)CC2=CC=CC=C2 NAME: '
               '2-[[6-methyl-4-oxo-5-(phenylmethyl)-1H-pyrimidin-2-yl]thio]hexanoic '
               'acid REASON: MISSED Contains ring structure(s), usually not '
               'characteristic of typical fatty acids\n'
               ' * SMILES: O1C(CCCC(O)C)=C(C(=C1CCC(O)=O)C(O)=O)C NAME: '
               '3-carboxy-4-methyl-5-(4-hydroxypentyl)-2-furanpropanoic acid '
               'REASON: MISSED Contains ring structure(s), usually not '
               'characteristic of typical fatty acids\n'
               ' * SMILES: '
               'C(CCCCCCCCC1C(CCCCCCCCCCC2C(CCCCCCCCCCCCCCCCC[C@@H](O)[C@H](C([O-])=O)CCCCCCCCCCCCCCCCCCCCCCCC)C2)C1)CCCCCCCCCC '
               'NAME: '
               '(2R)-2-[(1R)-1-hydroxy-18-{2-[10-(2-nonadecylcyclopropyl)decyl]cyclopropyl}octadecyl]hexacosanoate '
               'REASON: MISSED Carbon chain length 79 not in [4, 28]\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O(C(CCCCCCCCCCCC(O)=O)CCC)C(=O)CCCCCCCCCCCCCCC',
                                      'name': 'FAHFA(16:0/13-O-16:0)',
                                      'reason': 'Carbon chain length 32 for '
                                                'complexities beyond typical'},
                                  {   'smiles': 'O1C(CCCC)=CC=C1CCC2=CC(OC)=C(O)C=C2',
                                      'name': '2-Butyl-5-[2-(4-hydroxy-3-methoxyphenyl)ethyl]furan',
                                      'reason': 'No carboxylic acid group '
                                                'found'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OC(CCCCCCCC)CCCCCCCCC(O)=O',
                                      'name': '10-PAHSA',
                                      'reason': 'Carbon chain length 34 for '
                                                'complexities beyond typical'},
                                  {   'smiles': 'O=C(O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)C',
                                      'name': 'Fulvivirgamide B2',
                                      'reason': 'Carbon chain length 29 for '
                                                'complexities beyond typical'},
                                  {   'smiles': 'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCC(O)=O',
                                      'name': '(12Z,15Z,18Z,21Z,24Z)-triacontapentaenoic '
                                              'acid',
                                      'reason': 'Carbon chain length 30 for '
                                                'complexities beyond typical'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCC(=O)OC(CCCCC)CCCCCCCCCCCC(O)=O',
                                      'name': '13-PAHSA',
                                      'reason': 'Carbon chain length 34 for '
                                                'complexities beyond typical'},
                                  {   'smiles': 'OCCC(O)=O',
                                      'name': '3-hydroxypropionic acid',
                                      'reason': 'Carbon chain length 3 too '
                                                'short for typical fatty acid'},
                                  {   'smiles': 'CCC(C)CCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O',
                                      'name': '28-methyltriacontanoic acid',
                                      'reason': 'Carbon chain length 31 for '
                                                'complexities beyond typical'},
                                  {   'smiles': 'OC([C@@H](C[C@@H](C[C@@H](C[C@@H](C[C@@H](C[C@@H](C[C@@H](C[C@H](C)C(O)=O)C)C)C)C)C)C)C)CCCCCCCCCCCCCCCCC',
                                      'name': 'Hydroxyphthioceranic acid (C42)',
                                      'reason': 'Carbon chain length 42 for '
                                                'complexities beyond typical'},
                                  {   'smiles': 'O(C(CCCCCCCCCCCCCCC)=O)C(CCCCCCCCC)CCCCCCCC(=O)O',
                                      'name': '9-PAHSA',
                                      'reason': 'Carbon chain length 34 for '
                                                'complexities beyond typical'}],
    'attempt': 2,
    'reasoning': 'To improve the classification of fatty acids, we need to '
                 'address the false negative cases where the previous program '
                 "failed. Let's analyze the possible reasons for these "
                 'failures and think of improvements:\n'
                 '\n'
                 '1. **Ring Structure Detection**: The current code considers '
                 'any molecule with a ring structure as not a fatty acid. '
                 'However, some fatty acids can have cyclic substructures, '
                 'like Herbarin A, which contains a non-characteristic cyclic '
                 'acetal. We should not outright dismiss molecules with small '
                 'ring systems, especially if these are part of a larger '
                 'aliphatic chain. A solution might be to check if these rings '
                 'are small (5 or 6-membered) and not the predominant motif in '
                 'the molecule.\n'
                 '\n'
                 '2. **Carbon Chain Length**: The previous code was strict in '
                 'accepting only carbon chains between 4 and 28 atoms. '
                 'However, examples such as FAHFA(16:0/13-O-16:0) and some '
                 'complex derivatives included larger carbon chains, possibly '
                 'due to esterified or branched forms these may possess.\n'
                 '\n'
                 '3. **Multiple Substructural Features**: Some examples '
                 'contain multiple or complex functional groups. We need to '
                 'ensure a focus on an aliphatic backbone even when rings or '
                 'other groups are present. Perhaps evaluating connectivity or '
                 'predominant chain could help refine the logic.\n'
                 '\n'
                 '4. **Semantic Interpretation**: Clearing the ambiguity '
                 'around what truly constitutes a "fatty acid" beyond just '
                 'molecular features might require context-specific adjustment '
                 'in classifications.\n'
                 '\n'
                 'Essentially, our improvements will include more nuanced '
                 'handling of ring structures, tailored logic for chain '
                 'evaluation, and sensitivity to broader chemical contexts '
                 'when determining what pertains to the fatty acid class.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 2552,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9622926093514329,
    'f1': 0.9807840122982321,
    'accuracy': 0.9622926093514329,
    'negative_predictive_value': 0.0}