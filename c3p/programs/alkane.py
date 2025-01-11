"""
Classifies: CHEBI:18310 alkane
"""
from rdkit import Chem

def is_alkane(smiles: str):
    """
    Determines if a molecule is an alkane based on its SMILES string.
    An alkane is characterized by the formula CnH2n+2 and consists of
    saturated carbon atoms and hydrogen atoms, without any rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an alkane, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for presence of only carbon and hydrogen atoms
    elements = {atom.GetAtomicNum() for atom in mol.GetAtoms()}
    if elements.difference({6, 1}):  # Atomic number 6 is C, 1 is H
        return False, "Contains atoms other than carbon and hydrogen"
    
    # Check if the molecule is saturated (only single bonds)
    for bond in mol.GetBonds():
        if bond.GetBondTypeAsDouble() != 1.0:
            return False, "Contains unsaturated bonds (double or triple bonds present)"
    
    # Check for acyclic structure (cannot have rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Contains rings, not acyclic"
    
    # Correctly count carbons
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Correctly calculate the total hydrogen count, including implicit hydrogens
    hydrogen_count = sum(atom.GetTotalNumHs() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)

    # Calculate expected hydrogen count based on alkanes' CnH2n+2 rule
    expected_hydrogen_count = 2 * carbon_count + 2

    if hydrogen_count != expected_hydrogen_count:
        return False, f"Formula C{carbon_count}H{hydrogen_count} does not match CnH2n+2 (expected H{expected_hydrogen_count})"
    
    return True, "Molecule matches the definition of an alkane"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18310',
                          'name': 'alkane',
                          'definition': 'An acyclic branched or unbranched '
                                        'hydrocarbon having the general '
                                        'formula CnH2n+2, and therefore '
                                        'consisting entirely of hydrogen atoms '
                                        'and saturated carbon atoms.',
                          'parents': ['CHEBI:24632', 'CHEBI:33653'],
                          'xrefs': ['KEGG:C01371'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: CCCCCCCC(C)CCC NAME: 4-methylundecane '
               'REASON: MISSED Formula C12H0 does not match CnH2n+2 (expected '
               'H26)\n'
               ' * SMILES: CCC(C)CCC(C)CC(C)C(C)C NAME: '
               '2,3,5,8-tetramethyldecane REASON: MISSED Formula C14H0 does '
               'not match CnH2n+2 (expected H30)\n'
               ' * SMILES: C(CCCCCCCCCCCC(CCCCCCCCC)C)(CCCCCCCCCC)C NAME: '
               '10,22-Dimethyldotriacontane REASON: MISSED Formula C34H0 does '
               'not match CnH2n+2 (expected H70)\n'
               ' * SMILES: C(CCCCCCCCCCCCC)CCCCCCCCC(C)C NAME: '
               '2-Methyltetracosane REASON: MISSED Formula C25H0 does not '
               'match CnH2n+2 (expected H52)\n'
               ' * SMILES: C(CCCC(CC)C)(CCCC)C NAME: 3,7-dimethylundecane '
               'REASON: MISSED Formula C13H0 does not match CnH2n+2 (expected '
               'H28)\n'
               ' * SMILES: C(CC(C)C)(CCC)C NAME: 2,4-dimethylheptane REASON: '
               'MISSED Formula C9H0 does not match CnH2n+2 (expected H20)\n'
               ' * SMILES: CCC(CCCC(C)C)C(C)C NAME: 3-ethyl-2,7-dimethyloctane '
               'REASON: MISSED Formula C12H0 does not match CnH2n+2 (expected '
               'H26)\n'
               ' * SMILES: C(CC(C)C)(C(C)C)C NAME: 2,3,5-trimethylhexane '
               'REASON: MISSED Formula C9H0 does not match CnH2n+2 (expected '
               'H20)\n'
               ' * SMILES: CCC(C)CC(C)(C)CC NAME: 3,3,5-Trimethylheptane '
               'REASON: MISSED Formula C10H0 does not match CnH2n+2 (expected '
               'H22)\n'
               ' * SMILES: CCCCCCCCCCCCCCCC(C)CC NAME: 3-Methyloctadecane '
               'REASON: MISSED Formula C19H0 does not match CnH2n+2 (expected '
               'H40)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'Contains atoms other than '
                                               'carbon and hydrogen'}],
    'sample_false_negatives': [   {   'smiles': 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O',
                                      'name': 'octadecanedioic acid',
                                      'reason': 'Contains atoms other than '
                                                'carbon and hydrogen'}],
    'attempt': 5,
    'reasoning': 'The failure in the previous code was due to an incorrect '
                 'calculation of the hydrogen count, leading to the rejection '
                 'of valid alkane SMILES strings. The issue arose from the '
                 "malformed interpretation or handling of the molecule's "
                 'hydrogen atoms, especially considering implicit hydrogen '
                 "atoms within RDKit's molecule representation.\n"
                 '\n'
                 "Here's a systematic approach to resolving the identified "
                 'issues:\n'
                 '\n'
                 '1. **Fix Hydrogen Counting**: Ensure that all hydrogens, '
                 'including implicit ones, are counted correctly. Implicit '
                 'hydrogen atoms are common in RDKit, and we need to make sure '
                 'these are included in the count to accurately apply the '
                 'CnH2n+2 formula.\n'
                 '\n'
                 '2. **Structure Validation**: Confirm the structure is an '
                 'acyclic, saturated hydrocarbon by checking the bonds and '
                 'whether the molecule is not cyclic.\n'
                 '\n'
                 '3. **Testing**: Incorporate more robust testing or debug '
                 'print statements if necessary to confirm functionality when '
                 'handling example SMILES strings.\n'
                 '\n'
                 '### Updated Code',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 110,
    'num_false_positives': 3,
    'num_true_negatives': 142186,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.9734513274336283,
    'recall': 0.990990990990991,
    'f1': 0.9821428571428571,
    'accuracy': 0.9999718903724526,
    'negative_predictive_value': 0.9999929670082356}