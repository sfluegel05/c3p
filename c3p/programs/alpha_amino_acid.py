"""
Classifies: CHEBI:33704 alpha-amino acid
"""
from rdkit import Chem

def is_alpha_amino_acid(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid based on its SMILES string.
    An alpha-amino acid is defined as an amino acid in which the amino group 
    is located on the carbon atom at the position alpha to the carboxy group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxyl group pattern (C(=O)O)
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not carboxyl_matches:
        return False, "No carboxyl group found"

    # Look for an alpha-carbon bonded to both carboxyl (C) and amino group (N)
    for carboxyl_match in carboxyl_matches:
        carboxyl_c = carboxyl_match[0]  # C in 'C(=O)O'
        # Get neighbors of this C to find the alpha carbon
        neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(carboxyl_c).GetNeighbors() if nbr.GetAtomicNum() == 6]
        
        for alpha_c in neighbors:
            # Get neighbors of the alpha carbon and check if any is a nitrogen
            alpha_neighbors = [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(alpha_c).GetNeighbors()]
            if any(mol.GetAtomWithIdx(nbr).GetAtomicNum() == 7 for nbr in alpha_neighbors):  # N is atomic number 7
                return True, "Contains carboxyl group and amino group at alpha position"

    return False, "No amino group at alpha position found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33704',
                          'name': 'alpha-amino acid',
                          'definition': 'An amino acid in which the amino '
                                        'group is located on the carbon atom '
                                        'at the position alpha to the carboxy '
                                        'group.',
                          'parents': ['CHEBI:33709'],
                          'xrefs': ['KEGG:C00045', 'KEGG:C05167'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: CN[C@@H](C(C)C)C(O)=O NAME: '
               'N-methyl-L-valine REASON: MISSED No carboxyl group found\n'
               ' * SMILES: O=C(O)[C@@H](N)CC=1C=C(C(=CC1)O)OC NAME: '
               '3-O-methyldopa REASON: MISSED No carboxyl group found\n'
               ' * SMILES: CC(C)(CCO)SC[C@H](N)C(O)=O NAME: felinine REASON: '
               'MISSED No carboxyl group found\n'
               ' * SMILES: C[C@H](N)C(O)=O NAME: L-alanine REASON: MISSED No '
               'carboxyl group found\n'
               ' * SMILES: N[C@H](C(=O)O)CC=1C2=C(NC1)C=CC(=C2)O NAME: '
               '5-hydroxy-L-tryptophan REASON: MISSED No carboxyl group found\n'
               ' * SMILES: '
               '[H][C@](N)(CS[C@H](\\C=C\\C=C\\C=C/C\\C=C/CCC(O)=O)[C@@H](O)CCCC(O)=O)C(O)=O '
               'NAME: 18-carboxy-19,20-dinor-leukotriene E4 REASON: MISSED No '
               'carboxyl group found\n'
               ' * SMILES: C[C@H](O)[C@H](C)[C@H](N)C(O)=O NAME: '
               '(4S)-4-hydroxy-L-isoleucine REASON: MISSED No carboxyl group '
               'found\n'
               ' * SMILES: CC(N)(Cc1cnc[nH]1)C(O)=O NAME: '
               'alpha-methylhistidine REASON: MISSED No carboxyl group found\n'
               ' * SMILES: '
               'FC=1C(N2CCN(CC2)CC=3C=C([C@]4(N[C@@](CC5=C4NC=6C5=CC=CC6)(C(O)=O)[H])[H])C=CC3OC)=CC=CC1 '
               'NAME: trans-Ned 19 REASON: MISSED No carboxyl group found\n'
               ' * SMILES: N[C@@H](CC1=CC=C(F)C=C1)C(O)=O NAME: '
               '4-fluorophenyl-L-alanine REASON: MISSED No carboxyl group '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'No amino group at alpha '
                                               'position found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'No amino group at alpha '
                                               'position found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'No carboxyl group found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'No amino group at alpha '
                                               'position found'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'No amino group at alpha '
                                               'position found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'No carboxyl group found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'No carboxyl group found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No amino group at alpha '
                                               'position found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'No amino group at alpha '
                                               'position found'},
                                 {   'smiles': 'O1[C@H]([C@H](NC(=O)C)[C@@H](O)C[C@]1(O[C@@H]2[C@@H](O)[C@@H](O[C@@H]([C@@H]2O)CO)O[C@H]3[C@H](O)[C@@H](O)[C@@H](O[C@@H]3CO)O)C(O)=O)[C@H](O)[C@H](O[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](OC(=O)C)[C@H](O)CO)C(O)=O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-6-[(1S,2R)-2-[(2S,4S,5R,6R)-5-acetamido-6-[(1R,2R)-1-acetyloxy-2,3-dihydroxypropyl]-2-carboxy-4-hydroxyoxan-2-yl]oxy-1,3-dihydroxypropyl]-2-[(2R,3S,4S,5R,6S)-3,5-dihydroxy-2-(hydroxymethyl)-6-[(2R,3S,4R,5R,6R)-4,5,6-trihydroxy-2-(hydroxymethyl)oxan-3-yl]oxyoxan-4-yl]oxy-4-hydroxyoxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No amino group at alpha '
                                               'position found'}],
    'sample_false_negatives': [   {   'smiles': 'C(CC(=O)[O-])[C@@]1([H])C(=N[C@](CCC(=O)[O-])([H])C(=N1)O)O',
                                      'name': 'cyclo(Glu-Glu)',
                                      'reason': 'No amino group at alpha '
                                                'position found'},
                                  {   'smiles': 'O=C1[C@H]2[C@@](CC[C@H]3[C@H]2[C@@H](O)[C@H](O)[C@H]4[C@@H](O)[C@H](O)CC[C@]34C)(C)[C@H](C1)[C@H]([C@@H]5OC(O)C[C@@H](C5)C(C)C)C',
                                      'name': 'Contignasterol',
                                      'reason': 'No carboxyl group found'}],
    'attempt': 2,
    'reasoning': 'The previous attempt at classifying alpha-amino acids based '
                 'on their SMILES string failed because of misunderstanding '
                 'the carboxyl group representation in the SMILES format. The '
                 'SMILES for a neutral carboxyl group is typically "C(=O)O" '
                 'rather than "C(=O)[O-]". The latter represents a '
                 'deprotonated carboxylate ion, which is not applicable to all '
                 'alpha-amino acids.\n'
                 '\n'
                 'Additionally, the amino group detection might not correctly '
                 "identify the alpha carbon's relationship with the amino "
                 'group and carboxyl group due to inappropriate pattern '
                 'matching logic.\n'
                 '\n'
                 'To correctly identify an alpha-amino acid:\n'
                 '- The SMILES string must contain a carboxyl group "C(=O)O".\n'
                 '- There should be an amino group "N" attached to the carbon '
                 'alpha to the carboxyl carbon. The alpha carbon is the first '
                 'carbon attached to the carboxyl carbon.\n'
                 '\n'
                 'We can resolve these issues by:\n'
                 '- Correctly identifying the SMARTS pattern for the carboxyl '
                 'and amino groups.\n'
                 '- Mapping the alpha carbon appropriately and ensuring the '
                 'nitrogen is directly linked to this carbon via the SMILES '
                 'patterns.\n'
                 '\n'
                 "Here's the revised code with these changes incorporated:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 864,
    'num_false_positives': 100,
    'num_true_negatives': 840,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.8962655601659751,
    'recall': 0.9976905311778291,
    'f1': 0.9442622950819674,
    'accuracy': 0.9435215946843853,
    'negative_predictive_value': 0.997624703087886}