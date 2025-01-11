"""
Classifies: CHEBI:25413 monounsaturated fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import rdchem

def is_monounsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a monounsaturated fatty acid (MUFA) based on its SMILES string.
    MUFAs have one double or triple bond in the fatty acid chain and a carboxylic acid group.
    They must have a linear or cyclic configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monounsaturated fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define carboxylic acid pattern to confirm presence
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Define patterns for chain double and triple bonds that are NOT in rings
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")  # Represents C=C
    triple_bond_pattern = Chem.MolFromSmarts("[#6]#[#6]")  # Represents C#C

    # Detect number of double bonds and ensure they're not counted from carboxylic acid
    num_double_bonds = 0
    num_triple_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == rdchem.BondType.DOUBLE and not bond.IsInRing():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Ensure carbon backbone, ignore carboxylic C=O
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_double_bonds += 1
        elif bond.GetBondType() == rdchem.BondType.TRIPLE and not bond.IsInRing():
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            # Ensure carbon backbone
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                num_triple_bonds += 1

    # Count total unsaturations in chain
    num_unsaturations = num_double_bonds + num_triple_bonds

    # A MUFA should present exactly one instance of unsaturation in chain
    if num_unsaturations != 1:
        return False, f"Found {num_unsaturations} unsaturations in chain, need exactly one"
    
    return True, "Molecule is a monounsaturated fatty acid (one non-ring double or triple bond in the chain with a carboxylic group)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:25413',
                          'name': 'monounsaturated fatty acid',
                          'definition': 'Any fatty acid with one double or '
                                        'triple bond in the fatty acid chain '
                                        'and singly bonded carbon atoms in the '
                                        'rest of the chain. MUFAs have '
                                        'positive effects on the '
                                        'cardiovascular system, and in '
                                        'diabetes treatment.',
                          'parents': ['CHEBI:27208'],
                          'xrefs': ['PMID:10584045', 'PMID:12936956'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.038461538461538464 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: OC(=O)CCCCCCCCCCCC[C@H]1CCC=C1 NAME: '
               '(S)-chaulmoogric acid REASON: CORRECT Molecule is a '
               'monounsaturated fatty acid (one double or triple bond in the '
               'chain with carboxylic group)\n'
               ' * SMILES: CCCCCCCCC1=C(CCCCCCCC(O)=O)C1 NAME: sterculic acid '
               'REASON: CORRECT Molecule is a monounsaturated fatty acid (one '
               'double or triple bond in the chain with carboxylic group)\n'
               'False positives: NONE\n'
               'False negatives: SMILES: OCCCCCCCCCCCCC\\C=C\\C(O)=O NAME: '
               '(2E)-16-hydroxyhexadec-2-enoic acid REASON: MISSED Found 2 '
               'unsaturations in chain, need exactly one\n'
               ' * SMILES: [H]C(CCCCCCCC)=C([H])CCCCCCCC(O)=O NAME: '
               'octadec-9-enoic acid REASON: MISSED Found 2 unsaturations in '
               'chain, need exactly one\n'
               ' * SMILES: C(\\CCCCCCCC(=O)O)=C\\CCCCCCCC(=O)[H] NAME: '
               '18-oxooleic acid REASON: MISSED Found 3 unsaturations in '
               'chain, need exactly one\n'
               ' * SMILES: CCCCCCCC[C@H](O)\\C=C\\CCCCCCC(O)=O NAME: '
               '(8E,10S)-10-hydroxy-8-octadecenoic acid REASON: MISSED Found 2 '
               'unsaturations in chain, need exactly one\n'
               ' * SMILES: C(CCCCCCC/C=C\\C([C@H]1[C@H](CCCCC)O1)O)(=O)O NAME: '
               '11-hydroxy-(12S,13S)-epoxy-(9Z)-octadecenoic acid REASON: '
               'MISSED Found 2 unsaturations in chain, need exactly one\n'
               ' * SMILES: CCCC\\C=C/CCCCCCCCCC(O)=O NAME: '
               '(Z)-hexadec-11-enoic acid REASON: MISSED Found 2 unsaturations '
               'in chain, need exactly one\n'
               ' * SMILES: CCCCCC\\C=C/CCCC(O)=O NAME: cis-5-dodecenoic acid '
               'REASON: MISSED Found 2 unsaturations in chain, need exactly '
               'one\n'
               ' * SMILES: OCCCCCCCCCCCCCCCCC\\C=C\\C(O)=O NAME: '
               '(2E)-20-hydroxyicos-2-enoic acid REASON: MISSED Found 2 '
               'unsaturations in chain, need exactly one\n'
               ' * SMILES: C(CCCC(O)=O)CCCCCCCCCCCC#C NAME: octadec-17-ynoic '
               'acid REASON: MISSED Found 2 unsaturations in chain, need '
               'exactly one\n'
               ' * SMILES: [H]C(C)=CCCCCCCCCC(O)=O NAME: 10-dodecenoic acid '
               'REASON: MISSED Found 2 unsaturations in chain, need exactly '
               'one\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'Found 0 unsaturations in '
                                               'chain, need exactly one'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'Found 0 unsaturations in '
                                               'chain, need exactly one'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'Found 0 unsaturations in '
                                               'chain, need exactly one'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'Found 2 unsaturations in '
                                               'chain, need exactly one'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'Found 0 unsaturations in '
                                               'chain, need exactly one'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Found 0 unsaturations in '
                                               'chain, need exactly one'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'Found 0 unsaturations in '
                                               'chain, need exactly one'}],
    'sample_false_negatives': [   {   'smiles': 'OC(=O)CCCCCCCCCCCC[C@H]1CCC=C1',
                                      'name': '(S)-chaulmoogric acid',
                                      'reason': 'Found 0 unsaturations in '
                                                'chain, need exactly one'},
                                  {   'smiles': 'CCCCCCCCC1=C(CCCCCCCC(O)=O)C1',
                                      'name': 'sterculic acid',
                                      'reason': 'Found 0 unsaturations in '
                                                'chain, need exactly one'},
                                  {   'smiles': 'CCCCCCCCC1=C(CCCCCCC(O)=O)C1',
                                      'name': 'malvalic acid',
                                      'reason': 'Found 0 unsaturations in '
                                                'chain, need exactly one'},
                                  {   'smiles': 'OC(=O)CCCCCCCCCC[C@@H]1CCC=C1',
                                      'name': '(R)-hydnocarpic acid',
                                      'reason': 'Found 0 unsaturations in '
                                                'chain, need exactly one'},
                                  {   'smiles': 'OC(=O)CCCCCCCCCCCCC1CCC=C1',
                                      'name': 'chaulmoogric acid',
                                      'reason': 'Found 0 unsaturations in '
                                                'chain, need exactly one'},
                                  {   'smiles': 'O=C(CC(=O)O)CC/C=C/[C@@H](C/C=C\\CCCCC)O',
                                      'name': '(6E,8R,10Z)-8-hydroxy-3-oxohexadeca-6,10-dienoic '
                                              'acid',
                                      'reason': 'Found 2 unsaturations in '
                                                'chain, need exactly one'},
                                  {   'smiles': 'CCCCCCCC\\C=C/CCCCCCCC(=O)O[C@H](COC(=O)CCCC\\C=C/C\\C=C/C\\C=C/CCCCC)COP(O)(O)=O',
                                      'name': '1-(gamma-linolenoyl)-2-oleoyl-sn-glycero-3-phosphate',
                                      'reason': 'Found 4 unsaturations in '
                                                'chain, need exactly one'}],
    'attempt': 4,
    'reasoning': 'To address the issue with the previous code, several aspects '
                 'should be clarified and improved:\n'
                 '\n'
                 "1. **The Definition of 'Monounsaturated':** The previous "
                 'code mainly checks for the presence of a single non-ring '
                 'double or triple bond. However, it neglects to ensure that '
                 'the chain is linear or cyclic with one unsaturation that '
                 'fits the definition of a monounsaturated fatty acid.\n'
                 '\n'
                 '2. **Chain Consideration:** We need to ensure that the '
                 'double or triple bond is located in a long, primarily '
                 'carbon-based chain with a terminal carboxylic acid group. \n'
                 '\n'
                 '3. **Ignore Confusing Errors:** Ensure functional groups '
                 'that may masquerade as double bonds (such as carbonyl groups '
                 'in carboxylic acids) do not contribute to miscounting of '
                 'bond types.\n'
                 '\n'
                 '4. **Carboxylic Acid Check:** Verify a direct connection '
                 'between the unsaturated chain and the carboxylic acid '
                 'group.\n'
                 '\n'
                 "Let's incorporate the improvements and define the "
                 '`is_monounsaturated_fatty_acid` function to adequately '
                 'classify these molecules:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 198,
    'num_false_positives': 100,
    'num_true_negatives': 1597,
    'num_false_negatives': 7,
    'num_negatives': None,
    'precision': 0.6644295302013423,
    'recall': 0.9658536585365853,
    'f1': 0.7872763419483102,
    'accuracy': 0.9437434279705573,
    'negative_predictive_value': 0.9956359102244389}