"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
"""
Classifies: 2-hydroxy fatty acid
Definition: Any fatty acid with a hydroxy functional group in the alpha- or 2-position.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxy fatty acid based on its SMILES string.
    A 2-hydroxy fatty acid has a carboxylic acid group and a hydroxyl group at the 2-position,
    along with a sufficiently long carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find the carboxylic acid group
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxylic_acid_matches:
        return False, "No carboxylic acid group found"

    # Iterate over all carboxylic acid groups to find one with a hydroxyl at the 2-position
    for match in carboxylic_acid_matches:
        carboxylic_acid_atom = match[0]  # The carbon atom of the carboxylic acid group

        # Get the atom connected to the carboxylic acid carbon (the alpha carbon)
        alpha_carbon = None
        for neighbor in mol.GetAtomWithIdx(carboxylic_acid_atom).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon atom
                alpha_carbon = neighbor.GetIdx()
                break

        if alpha_carbon is None:
            continue  # No alpha carbon found, skip to the next match

        # Check if the alpha carbon has a hydroxyl group
        for neighbor in mol.GetAtomWithIdx(alpha_carbon).GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1:  # Hydroxyl group
                # Check for a sufficiently long carbon chain (fatty acid)
                carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
                if carbon_count >= 4:
                    return True, "Contains a carboxylic acid group, a hydroxyl group at the 2-position, and a sufficiently long carbon chain"

    return False, "No hydroxyl group at the 2-position relative to the carboxylic acid group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:10283',
                          'name': '2-hydroxy fatty acid',
                          'definition': 'Any fatty acid with a hydroxy '
                                        'functional group in the alpha- or '
                                        '2-position.',
                          'parents': ['CHEBI:24654', 'CHEBI:49302'],
                          'xrefs': [   'KEGG:C05102',
                                       'PMID:6423633',
                                       'PMID:9826358'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: CCCCCCCCCCCCCC(C)C(O)C(O)=O NAME: '
               '2-hydroxy-3-methylhexadecanoic acid REASON: MISSED No hydroxyl '
               'group at the 2-position relative to the carboxylic acid group\n'
               ' * SMILES: CCCCCCCCCC[C@H](O)C(O)=O NAME: (S)-2-hydroxylauric '
               'acid REASON: MISSED No hydroxyl group at the 2-position '
               'relative to the carboxylic acid group\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCC(O)C(O)=O NAME: '
               '(13Z,16Z,19Z,22Z)-2-hydroxyoctacosatetraenoic acid REASON: '
               'MISSED No hydroxyl group at the 2-position relative to the '
               'carboxylic acid group\n'
               ' * SMILES: CC(C)C[C@@H](O)C(O)=O NAME: '
               '(R)-2-hydroxy-4-methylpentanoic acid REASON: MISSED No '
               'hydroxyl group at the 2-position relative to the carboxylic '
               'acid group\n'
               ' * SMILES: CC(C)C[C@H](O)C(O)=O NAME: '
               '(S)-2-hydroxy-4-methylpentanoic acid REASON: MISSED No '
               'hydroxyl group at the 2-position relative to the carboxylic '
               'acid group\n'
               ' * SMILES: CCCCCCCCCCCCCC[C@@H](O)C(O)=O NAME: '
               '(R)-2-hydroxyhexadecanoic acid REASON: MISSED No hydroxyl '
               'group at the 2-position relative to the carboxylic acid group\n'
               ' * SMILES: CCCCCCCCCCCC(O)C(O)=O NAME: 2-hydroxytridecanoic '
               'acid REASON: MISSED No hydroxyl group at the 2-position '
               'relative to the carboxylic acid group\n'
               ' * SMILES: CCCCCCCCCCCCC(O)C(O)=O NAME: 2-hydroxymyristic acid '
               'REASON: MISSED No hydroxyl group at the 2-position relative to '
               'the carboxylic acid group\n'
               ' * SMILES: OC(C(CC)C)C(O)=O NAME: 2-hydroxy-3-methylpentanoic '
               'acid REASON: MISSED No hydroxyl group at the 2-position '
               'relative to the carboxylic acid group\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCCCCCC[C@@H](O)C(O)=O NAME: '
               '(R)-2-hydroxylignoceric acid REASON: MISSED No hydroxyl group '
               'at the 2-position relative to the carboxylic acid group\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)C[C@H](OC(=O)CCCCCCC/C=C\\CCCCCC)COC(=O)CCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC',
                                     'name': 'TG(22:1(13Z)/16:1(9Z)/18:4(6Z,9Z,12Z,15Z))',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O1C2C(C3(C(C4C(C5(C(CC4)CC(OC6OC(C(OC7OC(C(O)C(O)C7O)C)C(O)C6O)COC8OC(C(O)C(O)C8O)CO)CC5)C)CC3)C2)C)C(C19OCC(CC9)C)C',
                                     'name': 'Desglucoparillin',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': '[C@H]1(O)C(O)O[C@H](C=O)[C@H]([C@@H]1O)O',
                                     'name': '6-dehydro-D-glucose',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C1=CC(=CC=C1NC2=C(C=NC=C2)S(=O)(=O)N)Cl',
                                     'name': '4-(4-chloroanilino)-3-pyridinesulfonamide',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'C=1([NH+]=C(C=C(N1)NCC=2ON=C(N2)C(N)=O)C(C)C)N',
                                     'name': '2-amino-4-{[(3-carbamoyl-1,2,4-oxadiazol-5-yl)methyl]amino}-6-isopropylpyrimidin-1-ium',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'O([C@@H]1[C@@H](O[C@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)[C@H](O[C@@H]([C@@H]1O)CO)O)[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O)CO',
                                     'name': '(2R,3S,4S,5S,6R)-2-[(2S,3R,4S,5S,6R)-2,5-Dihydroxy-6-(hydroxymethyl)-3-[(2R,3R,4S,5S,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxyoxan-4-yl]oxy-6-(hydroxymethyl)oxane-3,4,5-triol',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)S(=O)(=O)NC2=CC3=C(C=C2)O[C@@H]4[C@H]3C[C@@H](O[C@@H]4CO)CC(=O)O',
                                     'name': '2-[(1R,3R,4aS,9aR)-1-(hydroxymethyl)-6-[(4-methoxyphenyl)sulfonylamino]-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]acetic '
                                             'acid',
                                     'reason': 'No hydroxyl group at the '
                                               '2-position relative to the '
                                               'carboxylic acid group'},
                                 {   'smiles': 'O=C(C=1NC=CC1)[C@@H]2[C@@H](C=C[C@H]3[C@H]2CCC3)C=CC=C([C@@H]4O[C@@H]([C@H](C(=O)O)C)CC[C@@H]4C)CC',
                                     'name': 'A83094A',
                                     'reason': 'No hydroxyl group at the '
                                               '2-position relative to the '
                                               'carboxylic acid group'},
                                 {   'smiles': 'C[C@H](O)[C@H](O)[C@@H](O)C(=O)COP(O)(O)=O',
                                     'name': 'L-rhamnulose 1-phosphate',
                                     'reason': 'No carboxylic acid group '
                                               'found'},
                                 {   'smiles': 'COc1ccc(\\C=C\\C(=O)N[C@H](CCCNC(N)=N)C(O)=O)cc1Br',
                                     'name': 'subereamine A',
                                     'reason': 'No hydroxyl group at the '
                                               '2-position relative to the '
                                               'carboxylic acid group'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'The previous code failed to correctly identify 2-hydroxy '
                 'fatty acids because the SMARTS pattern used to detect the '
                 'hydroxyl group at the 2-position relative to the carboxylic '
                 'acid group was too restrictive. Specifically, the pattern '
                 '`[CX3](=[OX1])[OX2H1][CX4H1,CX3H1][OX2H1]` assumes that the '
                 'hydroxyl group is directly attached to the carbon adjacent '
                 'to the carboxylic acid group, which is not always the case. '
                 'Additionally, the pattern does not account for the '
                 'possibility of branching or double bonds in the carbon '
                 'chain.\n'
                 '\n'
                 'To improve the program, we need to:\n'
                 '1. Use a more flexible SMARTS pattern that can detect the '
                 'hydroxyl group at the 2-position relative to the carboxylic '
                 'acid group, even if there are branches or double bonds in '
                 'the carbon chain.\n'
                 '2. Ensure that the hydroxyl group is correctly identified as '
                 'being at the 2-position relative to the carboxylic acid '
                 'group, regardless of the structure of the carbon chain.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 47,
    'num_false_positives': 100,
    'num_true_negatives': 14624,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3197278911564626,
    'recall': 1.0,
    'f1': 0.4845360824742268,
    'accuracy': 0.9932299776589263,
    'negative_predictive_value': 1.0}