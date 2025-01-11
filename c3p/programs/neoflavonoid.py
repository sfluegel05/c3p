"""
Classifies: CHEBI:71971 neoflavonoid
"""
from rdkit import Chem

def is_neoflavonoid(smiles: str):
    """
    Determines if a molecule is a neoflavonoid based on its SMILES string.
    A neoflavonoid is any 1-benzopyran with an aryl substituent at position 4.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a neoflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define more flexible SMARTS pattern for 1-benzopyran
    benzopyran_patterns = [
        Chem.MolFromSmarts("c1cc2occcc2c1"), # General 1-benzopyran
        Chem.MolFromSmarts("c1cc2[cH1]occc2c1"), # Allow for aromatic and hydrogen-bonding
        Chem.MolFromSmarts("C1=COc2ccccc2C1"), # Consider benzopyran with variation
    ]
    
    # Check for presence of 1-benzopyran structure
    for pattern in benzopyran_patterns:
        if mol.HasSubstructMatch(pattern):
            # Identify atoms that match the 1-benzopyran pattern
            benzopyran_matches = mol.GetSubstructMatches(pattern)

            # For each match, check if there's an aryl ring (aromatic six-membered ring) at position 4
            for match in benzopyran_matches:
                possible_aryl_positions = [4, 5]  # Possible indexing changes in different substructure mappings
                for pos_idx in possible_aryl_positions:
                    aryl_position_atom = mol.GetAtomWithIdx(match[pos_idx])
                    
                    # Scan the neighbors for an aryl group
                    for neighbor in aryl_position_atom.GetNeighbors():
                        # Check if this neighbor forms part of an aromatic six-membered ring
                        if neighbor.GetIsAromatic() and mol.GetRingInfo().IsAtomInRingOfSize(neighbor.GetIdx(), 6):
                            return True, "Aryl substituent found at or near position 4"
    
    return False, "No aryl substitution at or near position 4 in the 1-benzopyran"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:71971',
                          'name': 'neoflavonoid',
                          'definition': 'Any 1-benzopyran with an aryl '
                                        'substituent at position 4. The term '
                                        'was originally restricted to natural '
                                        'products, but is now also used to '
                                        'describe semi-synthetic and fully '
                                        'synthetic compounds.',
                          'parents': [   'CHEBI:26004',
                                         'CHEBI:38443',
                                         'CHEBI:72544'],
                          'xrefs': ['Wikipedia:Neoflavonoid'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'CC1=C(C=CC2=C1OC(=O)C=C2C3=CC=CC=C3)OCC(=O)O NAME: '
               '2-[(8-methyl-2-oxo-4-phenyl-1-benzopyran-7-yl)oxy]acetic acid '
               'REASON: MISSED No 1-benzopyran structure found\n'
               ' * SMILES: '
               'O=C(OC)C1=C(C(O)=CC(=C1)OC)C(=O)C2=C(O)C(C=3C4=C(OC(C3)(C)C)C(=CC(=C4)O)CC=C(C)C)=C(C)C=C2O '
               'NAME: Pestaloficiol L REASON: MISSED No 1-benzopyran structure '
               'found\n'
               ' * SMILES: '
               'O1[C@]2(O)[C@@]([C@H](C3=C1C(=C(O)C4=C3O[C@@H](CC4=O)C5=CC=CC=C5)C)C6=CC=CC=C6)(C(=O)C(C(=O)C2(C)C)(C)C)[H] '
               'NAME: Leucadenone B REASON: MISSED No 1-benzopyran structure '
               'found\n'
               ' * SMILES: '
               'O1C2=C(C(O)=C(C(O)=C2CC=C(C)C)C(=O)C(C)C)C(C3=CC=CC=C3)=CC1=O '
               'NAME: Mesuol REASON: MISSED No 1-benzopyran structure found\n'
               ' * SMILES: '
               'O1C(CC=2C1=C(C(O)=C3C2OC(=O)C=C3C4=CC=CC=C4)C(=O)C(C)C)C(O)(C)C '
               'NAME: '
               '1,2-Dihydro-5-hydroxy-2-(1-hydroxy-1-methylethyl)-4-(isobutyryl)-6-phenylfurano[2,3-h][1]benzopyran-8-one '
               'REASON: MISSED No 1-benzopyran structure found\n'
               ' * SMILES: O1C=2C(C(=CC1)C3=CC=CC=C3)=CC(O)=C(OC)C2OC NAME: '
               'Kuhlmannene REASON: MISSED No 1-benzopyran structure found\n'
               ' * SMILES: '
               'CC1=CC(N2C3=C1C=CC=C3C4(C2=O)C5=C(C=C(C=C5)O)OC(=C4C#N)N)(C)C '
               'NAME: LSM-26505 REASON: MISSED No 1-benzopyran structure '
               'found\n'
               ' * SMILES: '
               'O1C(C=CC=2C1=C(C(O)=C3C2OC(=O)C=C3C4=CC=CC=C4)C(=O)C(C)C)(CO)C '
               'NAME: '
               '5-hydroxy-8-(hydroxymethyl)-8-methyl-6-(2-methylpropanoyl)-4-phenyl-2H,8H-pyrano[2,3-f]chromen-2-one '
               'REASON: MISSED No 1-benzopyran structure found\n'
               ' * SMILES: O1C=2C(C(CC1=O)C3=CC=C(OC)C=C3)=C(O)C=C(O)C2 NAME: '
               '4-(4-Methoxyphenyl)-5,7-dihydroxychroman-2-one REASON: MISSED '
               'No 1-benzopyran structure found\n'
               ' * SMILES: '
               'O=C(C1=C(O)C(=CC=C1O)[C@@H]2C3=C(O[C@@H](C2)C)C=CC=C3OC)CCC '
               'NAME: Nodulisporin H REASON: MISSED No 1-benzopyran structure '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O([C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)CO)CO)[C@H]3[C@@H](O)[C@H](OC(O)[C@@H]3O)CO',
                                     'name': 'beta-D-Glcp-(1->2)-beta-D-Glcp-(1->3)-D-Galp',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': 'CCS(=O)(=O)NCC[C@@H]1CC[C@@H]([C@H](O1)CO)NC(=O)NC2=CC(=CC(=C2)Cl)Cl',
                                     'name': '1-(3,5-dichlorophenyl)-3-[(2S,3S,6S)-6-[2-(ethylsulfonylamino)ethyl]-2-(hydroxymethyl)-3-oxanyl]urea',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': 'C(C(O)=O)C/C=C\\C/C=C\\C\\C=C/C=C/C=C/[C@H]1[C@H](C/C=C\\CC)O1',
                                     'name': '(16S,17S)-epoxy-(4Z,7Z,10Z,12E,14E,19Z)-docosahexaenoic '
                                             'acid',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': 'OC[C@H]1O[C@H](O[C@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](CO)O[C@@H](O)[C@H](O)[C@H]3O)[C@H](O)[C@H]2O)[C@H](O)[C@@H](O)[C@H]1O',
                                     'name': 'alpha-D-Galp-(1->4)-beta-D-Galp-(1->4)-beta-D-Glcp',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': '[Li+].[Br-]',
                                     'name': 'lithium bromide',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': 'C=1(OC)C2=C(C=C3C[C@H]([C@](CC=4C=C(OC)C(OC)=C(C4C13)OC)(C)O)C)OCO2',
                                     'name': 'Besigomsin',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C=CC(=O)C2=CC=CN2)Cl',
                                     'name': '3-(2-chlorophenyl)-1-(1H-pyrrol-2-yl)-2-propen-1-one',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': 'O=C1OC(O)C(=C1C(O)C(C)C)C',
                                     'name': '5-Hydroxy-3-(1-hydroxy-2-methylpropyl)-4-methyl-2(5H)-furanone',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': '[H][C@@]1(COC(C)(C)[C@@]1([H])OC(=O)\\C=C\\c1ccccc1)c1c(O)ccc2C(=O)C[C@H](Oc12)c1ccccc1',
                                     'name': '(+)-tephrorin B',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)CC4=CC=CC(=C4)C(=O)O)C5=CC=CC=C5N2C)[C@@H](C)CO',
                                     'name': 'LSM-9341',
                                     'reason': 'No aryl substitution at or '
                                               'near position 4 in the '
                                               '1-benzopyran'}],
    'sample_false_negatives': [   {   'smiles': 'CC1=C(C=CC2=C1OC(=O)C=C2C3=CC=CC=C3)OCC(=O)O',
                                      'name': '2-[(8-methyl-2-oxo-4-phenyl-1-benzopyran-7-yl)oxy]acetic '
                                              'acid',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'O=C(OC)C1=C(C(O)=CC(=C1)OC)C(=O)C2=C(O)C(C=3C4=C(OC(C3)(C)C)C(=CC(=C4)O)CC=C(C)C)=C(C)C=C2O',
                                      'name': 'Pestaloficiol L',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'O1[C@]2(O)[C@@]([C@H](C3=C1C(=C(O)C4=C3O[C@@H](CC4=O)C5=CC=CC=C5)C)C6=CC=CC=C6)(C(=O)C(C(=O)C2(C)C)(C)C)[H]',
                                      'name': 'Leucadenone B',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'O1C2=C(C(O)=C(C(O)=C2CC=C(C)C)C(=O)C(C)C)C(C3=CC=CC=C3)=CC1=O',
                                      'name': 'Mesuol',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'O1C(CC=2C1=C(C(O)=C3C2OC(=O)C=C3C4=CC=CC=C4)C(=O)C(C)C)C(O)(C)C',
                                      'name': '1,2-Dihydro-5-hydroxy-2-(1-hydroxy-1-methylethyl)-4-(isobutyryl)-6-phenylfurano[2,3-h][1]benzopyran-8-one',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'O1C=2C(C(=CC1)C3=CC=CC=C3)=CC(O)=C(OC)C2OC',
                                      'name': 'Kuhlmannene',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'O1C(C=CC=2C1=C(C(O)=C3C2OC(=O)C=C3C4=CC=CC=C4)C(=O)C(C)C)(CO)C',
                                      'name': '5-hydroxy-8-(hydroxymethyl)-8-methyl-6-(2-methylpropanoyl)-4-phenyl-2H,8H-pyrano[2,3-f]chromen-2-one',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'O1C=2C(C(CC1=O)C3=CC=C(OC)C=C3)=C(O)C=C(O)C2',
                                      'name': '4-(4-Methoxyphenyl)-5,7-dihydroxychroman-2-one',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'O=C(C1=C(O)C(=CC=C1O)[C@@H]2C3=C(O[C@@H](C2)C)C=CC=C3OC)CCC',
                                      'name': 'Nodulisporin H',
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'},
                                  {   'smiles': 'CC1=CC(=CC2=C1C(CC3(O2)CC(NC(=S)N3)(C)C)C4=CC=C(C=C4)OC)O',
                                      'name': "7'-hydroxy-4'-(4-methoxyphenyl)-4,4,5'-trimethyl-2-spiro[1,3-diazinane-6,2'-3,4-dihydro-2H-1-benzopyran]thione",
                                      'reason': 'No aryl substitution at or '
                                                'near position 4 in the '
                                                '1-benzopyran'}],
    'attempt': 4,
    'reasoning': 'The previous program aimed to classify a SMILES string as a '
                 'neoflavonoid by checking two main criteria: the presence of '
                 'a 1-benzopyran structure and an aryl group attached at the '
                 '4-position. However, the program failed across various '
                 "instances, as indicated by the outcomes. Let's analyze the "
                 'challenges and possible improvements:\n'
                 '\n'
                 '### Analysis of Failures:\n'
                 '\n'
                 '1. **Pattern Definition for 1-benzopyran:**\n'
                 '   - The SMARTS pattern for the 1-benzopyran structure was '
                 'too specific or limited. It likely failed to capture the '
                 'structural diversity found in neoflavonoids like the '
                 'presence of heteroatoms and variations in aromaticity.\n'
                 '\n'
                 '2. **Substructure Matching Issues:**\n'
                 '   - It seems that the substructure matching technique is '
                 'not capturing all the molecules, especially with more '
                 'complex structures or fused systems.\n'
                 '\n'
                 '3. **Aryl Group Detection:**\n'
                 '   - The program checks for an aryl ring at the 4-position '
                 'based on indices in a SMARTS pattern. This rigid approach '
                 'may not account for all ways an aryl group might connect in '
                 'varied synthetics or natural products.\n'
                 '\n'
                 '### Improvement Plan:\n'
                 '\n'
                 '1. **Refine SMARTS Patterns:**\n'
                 '   - Use more flexible SMARTS patterns that can capture '
                 '1-benzopyran presence more reliably, considering alternative '
                 'aromatic ring systems and potential variations due to '
                 'substituted or modified moieties.\n'
                 '\n'
                 '2. **General Smiles Parsing:**\n'
                 '   - For highly variable molecules, complement SMARTS-based '
                 'substructure searches with more flexible graph or '
                 'scaffold-based methods. Additionally, manual inspection of '
                 'chemical rules could be beneficial.\n'
                 '\n'
                 '3. **Enhanced Positioning Logic:**\n'
                 '   - Beyond simple index checking, explore neighborhood '
                 'connectivity more robustly to verify the attachment of aryl '
                 'groups in flexible and variable bonding patterns.\n'
                 '\n'
                 '4. **Iterate on Specific Examples:**\n'
                 '   - Adjust the program based on how each specific failed '
                 'example deviates from expectations.\n'
                 '\n'
                 "Here's the revised program incorporating these enhancements:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 2,
    'num_false_positives': 23,
    'num_true_negatives': 142221,
    'num_false_negatives': 54,
    'num_negatives': None,
    'precision': 0.08,
    'recall': 0.03571428571428571,
    'f1': 0.04938271604938271,
    'accuracy': 0.9994588896697119,
    'negative_predictive_value': 0.9996204533473906}