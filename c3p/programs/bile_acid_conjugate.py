"""
Classifies: CHEBI:36249 bile acid conjugate
"""
"""
Classifies: Bile Acid Conjugate
"""

from rdkit import Chem

def is_bile_acid_conjugate(smiles: str):
    """
    Determines if a molecule is a bile acid conjugate based on its SMILES string.
    A bile acid conjugate is a bile acid attached to a hydrophilic group such as glycine, taurine,
    other amino acids, sulfate, glucuronic acid, glucose, other sugars, or coenzyme A.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid conjugate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for steroid nucleus: fused ring system of four rings (6-6-6-5)
    ri = mol.GetRingInfo()
    sssr = ri.AtomRings()
    if len(sssr) < 4:
        return False, "No steroid nucleus found (not enough rings)"

    # Get ring sizes
    ring_sizes = [len(ring) for ring in sssr]

    # Check for ring sizes 6-6-6-5
    six_membered_rings = ring_sizes.count(6)
    five_membered_rings = ring_sizes.count(5)
    if not (six_membered_rings >= 3 and five_membered_rings >= 1):
        return False, "No steroid nucleus found (incorrect ring sizes)"

    # Now check for hydrophilic conjugation groups
    # Define SMARTS patterns for conjugation groups

    # Glycine conjugation (amide bond to glycine)
    glycine_conj_smarts = '[NX3][CX3](=O)[CX4][CX3](=O)[O-,O]'
    glycine_conj = Chem.MolFromSmarts(glycine_conj_smarts)

    # Taurine conjugation (amide bond to taurine)
    taurine_conj_smarts = '[NX3][CX4][CX4][SX4](=O)(=O)[O-,O]'
    taurine_conj = Chem.MolFromSmarts(taurine_conj_smarts)

    # Sulfate conjugation (sulfate ester)
    sulfate_conj_smarts = '[OX2][SX4](=O)(=O)[O-,O]'
    sulfate_conj = Chem.MolFromSmarts(sulfate_conj_smarts)

    # Glucuronic acid conjugation (ester or glycosidic bond to glucuronic acid)
    glucuronic_acid_conj_smarts = 'O[C@H]1[C@@H](O)[C@H](O)[C@H](O)[C@@H]1C(=O)[O-,O]'
    glucuronic_acid_conj = Chem.MolFromSmarts(glucuronic_acid_conj_smarts)

    # Glucose conjugation (glycosidic bond)
    glucose_conj_smarts = 'O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@H]1O'
    glucose_conj = Chem.MolFromSmarts(glucose_conj_smarts)

    # General amino acid conjugation (amide bond to amino acid side chain)
    amino_acid_conj_smarts = '[NX3][CX3](=O)[CX4][CX3](=O)[O-,O]'
    amino_acid_conj = Chem.MolFromSmarts(amino_acid_conj_smarts)

    # List of conjugation patterns
    conjugation_patterns = [
        ('glycine', glycine_conj),
        ('taurine', taurine_conj),
        ('sulfate', sulfate_conj),
        ('glucuronic acid', glucuronic_acid_conj),
        ('glucose', glucose_conj),
        ('amino acid', amino_acid_conj)
    ]

    # Check for conjugation
    for name, pattern in conjugation_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, f"Contains {name} conjugation"

    # Also check for coenzyme A conjugation (simplified pattern)
    coa_conj_smarts = 'NC(=O)CCNC(=O)[CX4][CX3](=O)[O-,O]'
    coa_conj = Chem.MolFromSmarts(coa_conj_smarts)
    if mol.HasSubstructMatch(coa_conj):
        return True, "Contains coenzyme A conjugation"

    return False, "No hydrophilic conjugation found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36249',
                          'name': 'bile acid conjugate',
                          'definition': 'Any bile acid conjugated to a '
                                        'functional group that gives '
                                        'additional hydrophilicity or charge '
                                        'to the molecule. Molecules used for '
                                        'conjugation are: glycine, taurine '
                                        '(and other amino acids); sulfuric '
                                        "acid (for which the term ''sulfate'' "
                                        'may be used); glucuronic acid (for '
                                        "which the term ''glucuronate'' may be "
                                        'used); glucose and other uncharged '
                                        'sugars; and coenzyme A.',
                          'parents': ['CHEBI:36078'],
                          'xrefs': [   'PMID:19619630',
                                       'PMID:35224661',
                                       'PMID:36396870',
                                       'PMID:37071431'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: invalid syntax (<string>, line 1)\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'No steroid nucleus found (not '
                                               'enough rings)'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'No steroid nucleus found (not '
                                               'enough rings)'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'No hydrophilic conjugation '
                                               'found'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'No steroid nucleus found '
                                               '(incorrect ring sizes)'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'No steroid nucleus found (not '
                                               'enough rings)'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No steroid nucleus found (not '
                                               'enough rings)'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'No steroid nucleus found (not '
                                               'enough rings)'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'No steroid nucleus found (not '
                                               'enough rings)'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'No steroid nucleus found (not '
                                               'enough rings)'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'No steroid nucleus found '
                                               '(incorrect ring sizes)'}],
    'sample_false_negatives': [   {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3([C@H](C[C@@]2(C[C@@H](C1)O)[H])O)[H])(CC[C@]4([H])[C@@H](CCC(NCC(O)=O)=O)C)[H])C)[H])C',
                                      'name': 'glycoursodeoxycholic acid',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': 'OC1C2C3C(C(CC3)C(CCC(=O)NCC(O)=O)C)(CCC2C4(C(C1O)CC(O)CC4)C)C',
                                      'name': 'Glycine-|A-muricholic Acid',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': 'O[C@H]1[C@]2([C@]3([C@@]([C@](CC3)([C@@H](CCC(=O)NC(CO)C(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(C[C@H](O)CC4)[H])C)[H])C)[H])[H]',
                                      'name': 'Serine conjugated '
                                              'chenodeoxycholic acid',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': 'O[C@H]1C[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@](CC4)([C@@H](CCC(=O)N(CC(O)=O)C)C)[H])(CC3)C)[H])(CC2)[H])[H])(CC1)C)[H]',
                                      'name': 'N-((4R)-4-((3R,5S,9S,10S,13R,14S,17R)-3-Hydroxy-10,13-dimethylhexadecahydro-1H-cyclopenta[a]phenanthren-17-yl)pentanoyl)-N-methylglycine',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': 'OC1C2(C(C3C(C4(C(CC3O)CC(O)CC4)C)C1)CCC2C(CCC(=O)NC(CCCNC(=O)N)C(O)=O)C)C',
                                      'name': 'Cholylcitrulline',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': '[H][C@]12CC[C@@]3([H])[C@]4([H])CC[C@]([H])([C@H](C)CCC(=O)NCC(O)=O)[C@@]4(C)CC[C@]3([H])[C@@]1(C)CC[C@@H](O)C2',
                                      'name': 'glycolithocholic acid',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': 'OC1C2(C(C3C(C4(C(CC3O)CC(O)CC4)C)C1)CCC2C(CCC(=O)NC(C)C(O)=O)C)C',
                                      'name': 'Cholylalanine',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': 'O[C@H]1[C@]2([C@]3([C@@]([C@](CC3)([C@@H](CCC(=O)NC(CCC(O)=O)C(O)=O)C)[H])(CC[C@@]2([C@@]4([C@](C1)(C[C@H](O)CC4)[H])C)[H])C)[H])[H]',
                                      'name': 'Glutamate conjugated '
                                              'chenodeoxycholic acid',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': 'OC1C2(C(C3C(C4(C(CC3)CC(O)CC4)C)C1)CCC2C(CCC(=O)NC(CCC(=O)N)C(O)=O)C)C',
                                      'name': 'Deoxycholylglutamine',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'},
                                  {   'smiles': 'OC1C2(C(C3C(C4(C(CC3O)CC(O)CC4)C)C1)CCC2C(CCC(=O)NC(CCCCN)C(O)=O)C)C',
                                      'name': 'Lysine conjugated cholic acid',
                                      'reason': 'No hydrophilic conjugation '
                                                'found'}],
    'attempt': 2,
    'reasoning': 'The previous code failed due to a syntax error caused by an '
                 'unterminated triple-quoted string at line 123. This '
                 'indicates that the code had an incomplete multi-line string, '
                 'likely in a docstring or comment, which led to invalid '
                 'syntax. Additionally, the F1 score was 0, suggesting that '
                 'the function did not correctly identify any true positives. '
                 'This could be due to incorrect or insufficient substructure '
                 'patterns used to detect bile acid conjugates.\n'
                 '\n'
                 'To improve, we need to ensure that all triple-quoted strings '
                 'are properly closed to prevent syntax errors. We should also '
                 'enhance the detection of the bile acid core and conjugation '
                 'groups by carefully defining the SMARTS patterns that '
                 'accurately reflect the structures of bile acids and their '
                 'conjugates.\n'
                 '\n'
                 "In the new code, we'll:\n"
                 '\n'
                 '- Correct any syntax errors by ensuring all strings are '
                 'properly closed.\n'
                 '- Define a robust method to identify the steroid nucleus '
                 'characteristic of bile acids.\n'
                 '- Use comprehensive SMARTS patterns for the various '
                 'conjugation groups (glycine, taurine, sulfate, sugars, amino '
                 'acids).\n'
                 '- Ensure that the function correctly identifies the examples '
                 'provided.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 39,
    'num_false_positives': 100,
    'num_true_negatives': 53463,
    'num_false_negatives': 58,
    'num_negatives': None,
    'precision': 0.2805755395683453,
    'recall': 0.4020618556701031,
    'f1': 0.3305084745762712,
    'accuracy': 0.9970555348490495,
    'negative_predictive_value': 0.9989163132228471}