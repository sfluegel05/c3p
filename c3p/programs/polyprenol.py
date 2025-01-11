"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a prenol with more than one isoprene unit and a terminal hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for repeating isoprene units
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C-C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    if len(isoprene_matches) <= 1:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than 1"

    # Check for terminal hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing terminal hydroxyl group"

    # Check the position of the hydroxyl group to ensure it is at the end
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    terminal_oh = False
    for match in hydroxyl_matches:
        atom = mol.GetAtomWithIdx(match[0])
        neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
        # Ensure that OH is attached to a carbon (often the terminal position in polyprenols)
        if 6 in neighbors and len(neighbors) == 1:
            terminal_oh = True
            break

    if not terminal_oh:
        return False, "Hydroxyl group is not terminal"

    return True, "Contains more than one isoprene unit with a terminal hydroxyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26199',
                          'name': 'polyprenol',
                          'definition': 'Any member of the class of  prenols '
                                        'possessing the general formula '
                                        'H-[CH2C(Me)=CHCH2]nOH in which the '
                                        'carbon skeleton is composed of more '
                                        'than one isoprene units.',
                          'parents': ['CHEBI:26244'],
                          'xrefs': ['KEGG:C06081'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'Found 1 isoprene units, need '
                                               'more than 1'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'Found 0 isoprene units, need '
                                               'more than 1'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'Found 0 isoprene units, need '
                                               'more than 1'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'Found 0 isoprene units, need '
                                               'more than 1'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'Found 0 isoprene units, need '
                                               'more than 1'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'Missing terminal hydroxyl '
                                               'group'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'Found 0 isoprene units, need '
                                               'more than 1'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Found 0 isoprene units, need '
                                               'more than 1'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'Found 0 isoprene units, need '
                                               'more than 1'},
                                 {   'smiles': 'O1[C@H]([C@H](NC(=O)C)[C@@H](O)C[C@]1(O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O)C(O)=O)[C@H](O)[C@H](O[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](OC(=O)C)[C@H](O)CO)C(O)=O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-6-[(1S,2R)-2-[(2S,4S,5R,6R)-5-acetamido-6-[(1R,2R)-1-acetyloxy-2,3-dihydroxypropyl]-2-carboxy-4-hydroxyoxan-2-yl]oxy-1,3-dihydroxypropyl]-2-[(2R,3S,4S,5R,6S)-3,5-dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-4-hydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Found 0 isoprene units, need '
                                               'more than 1'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'To classify a molecule as belonging to the class of '
                 'polyprenols, we can create a function that checks for the '
                 "presence of certain structural features in the molecule's "
                 'SMILES string. A polyprenol is a type of prenol with the '
                 'general formula H-[CH2C(Me)=CHCH2]nOH, meaning it consists '
                 'of multiple isoprene units (C5H8) in which "n" must be more '
                 'than 1. Here are key characteristics to look for in a '
                 'polyprenol:\n'
                 '\n'
                 '1. A repeating isoprene unit structure, which in SMILES '
                 'could be represented as:\n'
                 '\n'
                 '   - C-C(=C)-C-C\n'
                 '   \n'
                 '2. A hydroxyl group (-OH) at the end of the chain.\n'
                 '\n'
                 '3. The molecule must contain more than one isoprene unit, '
                 'meaning the chain is longer than just a single isoprene.\n'
                 '\n'
                 'We can use RDKit, a cheminformatics library, to parse SMILES '
                 'strings, and then search for these patterns in the molecule. '
                 "Here's the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 27,
    'num_false_positives': 100,
    'num_true_negatives': 527,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2125984251968504,
    'recall': 1.0,
    'f1': 0.35064935064935066,
    'accuracy': 0.8470948012232415,
    'negative_predictive_value': 1.0}