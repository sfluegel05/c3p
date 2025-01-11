"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is characterized by a hydroxyl (-OH) group attached
    to a carbon which is attached to at least two hydrogen atoms, or a
    carbon with one other carbon and two hydrogen atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define an improved SMARTS pattern that captures a broader definition of primary alcohol
    # The main features of a primary alcohol: OH group not linked to carbon with more than one non-hydrogen substituent
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4;H2,H1][OH]")

    # Check the pattern match
    if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "SMILES indicates a primary alcohol structure"
    else:
        return False, "No primary alcohol structure found in SMILES"

# Example of usage:
# smiles = '...'  # Replace with actual SMILES string
# is_primary_alcohol(smiles)


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15734',
                          'name': 'primary alcohol',
                          'definition': 'A primary alcohol is a compound in '
                                        'which a hydroxy group, -OH, is '
                                        'attached to a saturated carbon atom '
                                        'which has either three hydrogen atoms '
                                        'attached to it or only one other '
                                        'carbon atom and two hydrogen atoms '
                                        'attached to it.',
                          'parents': ['CHEBI:30879'],
                          'xrefs': ['KEGG:C00226'],
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
               'COc1ccc([C@H](O)[C@H](CO)Oc2c(OC)cc(cc2OC)[C@@H]2OC[C@@H]3[C@H]2CO[C@H]3c2cc(OC)c(O)c(OC)c2)c(OC)c1O '
               'NAME: '
               "(-)-(7R,7'R,7''S,8S,8'S,8''S)-4',4''-dihydroxy-3,3',3'',5,5',5''-hexamethoxy-7,9':7',9-diepoxy-4,8''-oxy-8,8'-sesquineolignan-7'',9''-diol "
               'REASON: MISSED No primary alcohol structure found in SMILES\n'
               ' * SMILES: [H][C@@]1(CCC(C)=CC1)C(C)(C)O NAME: '
               '(S)-(-)-alpha-terpineol REASON: MISSED No primary alcohol '
               'structure found in SMILES\n'
               ' * SMILES: '
               'COc1cc(ccc1O)C(=O)CC(O)c1cc(\\C=C\\CO)cc(OC)c1O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: bidenlignaside A REASON: MISSED No primary alcohol '
               'structure found in SMILES\n'
               ' * SMILES: C\\C(CO)=C\\CC[C@]1(C)[C@H]2CC[C@H](C2)C1=C NAME: '
               'beta-santalol REASON: MISSED No primary alcohol structure '
               'found in SMILES\n'
               ' * SMILES: CCCCCCCCCCCCCCCCCC(=O)NCCO NAME: '
               'N-(octadecanoyl)ethanolamine REASON: MISSED No primary alcohol '
               'structure found in SMILES\n'
               ' * SMILES: '
               'OCCN1CCN(CC\\C=C2\\C3=C(SC4=C2C=C(C=C4)C(F)(F)F)C=CC=C3)CC1 '
               'NAME: cis-flupenthixol REASON: MISSED No primary alcohol '
               'structure found in SMILES\n'
               ' * SMILES: OCCCCCCCCCCCCC\\C=C\\C(O)=O NAME: '
               '(2E)-16-hydroxyhexadec-2-enoic acid REASON: MISSED No primary '
               'alcohol structure found in SMILES\n'
               ' * SMILES: CC(O)CCCCCCCCO NAME: 1,9-decanediol REASON: MISSED '
               'No primary alcohol structure found in SMILES\n'
               ' * SMILES: '
               '[C@@]1(C)(CO)CC[C@@H](C(=C1\\C=C\\C(=C\\C=C\\C(=C\\C(=O)O)\\C)\\C)C)O '
               'NAME: (4S)-4,16-dihydroxyretinoic acid REASON: MISSED No '
               'primary alcohol structure found in SMILES\n'
               ' * SMILES: '
               'COc1ccccc1Oc1c(NS(=O)(=O)c2ccc(cc2)C(C)(C)C)nc(nc1OCCO)-c1ncccn1 '
               'NAME: bosentan REASON: MISSED No primary alcohol structure '
               'found in SMILES\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': '[H][C@]12CCC[C@H](C)[C@@]1(C)C[C@@H](CC2)C(C)C',
                                     'name': 'eremophilane',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1NC=NC1)C(O)=O)[C@@H](NC(=O)[C@@H](N)CC2=CC=CC=C2)CCC(=O)N',
                                     'name': 'Phe-Gln-His',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'S(=O)(=O)(O[C@@H](C(=O)N[C@@H]1C(=O)N[C@H](C(=O)N[C@@H]2C(=O)N([C@@H]([C@H](CC)C)C(N([C@H](C(N[C@H](C(O[C@@H]1C)=O)[C@H](CC)C)=O)CC3=CC=CC=C3)C)=O)[C@H](OC)CC2)CCCN=C(N)N)COS(=O)(=O)O)O',
                                     'name': 'Micropeptin MZ1019',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'CC1=CC=C(C=C1)CS(=O)(=O)C2=NC=C(C(=N2)C(=O)NC3=NN=C(S3)C(C)C)Cl',
                                     'name': '5-chloro-2-[(4-methylphenyl)methylsulfonyl]-N-(5-propan-2-yl-1,3,4-thiadiazol-2-yl)-4-pyrimidinecarboxamide',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'P(OCC[N+](C)(C)C)(OC[C@H](OC(=O)CCCCCCC/C=C\\CCCCC)COC(=O)CCCCCCC/C=C\\CCCCC)([O-])=O',
                                     'name': 'PC(15:1(9Z)/15:1(9Z))',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CC(O)=O)C)[C@H](CC)C',
                                     'name': 'Asp-Ala-Ile',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'O=C1OC(=CC(=C1)O)[C@H]2C3=C(C(O)=CC=C3)C(=O)C[C@]2(O)C4=CC=CC=C4',
                                     'name': 'Wailupemycin D',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C3=CC=CC=C3O2)NC(=O)CSC4=NC5=CC=CC=C5N4',
                                     'name': '2-(1H-benzimidazol-2-ylthio)-N-(2-methoxy-3-dibenzofuranyl)acetamide',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'CC(=O)O[C@@H]1C[C@H]2C(C)(C)C(=O)C=C[C@]2(C)[C@H]2CC[C@@]3(C)[C@@H](C(=O)C=C3[C@]12C)c1ccoc1',
                                     'name': 'azadiradione',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'},
                                 {   'smiles': 'COC1=C(C=C(C=C1)N2C(=O)C3=CC=CC=C3N(C2=O)CC(=O)OCC4=CC=CC=C4)OC',
                                     'name': '2-[3-(3,4-dimethoxyphenyl)-2,4-dioxo-1-quinazolinyl]acetic '
                                             'acid (phenylmethyl) ester',
                                     'reason': 'No primary alcohol structure '
                                               'found in SMILES'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]1(CCC(C)=CC1)C(C)(C)O',
                                      'name': '(S)-(-)-alpha-terpineol',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'},
                                  {   'smiles': 'CO',
                                      'name': 'methanol',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'},
                                  {   'smiles': 'C[C@@](O)(CCOP(O)(O)=O)CC(O)=O',
                                      'name': '(R)-5-phosphomevalonic acid',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'},
                                  {   'smiles': 'O(C([H])([H])[H])[H]',
                                      'name': 'methanol-d1',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'},
                                  {   'smiles': '[H][C@]1(CCC(C)=CC1)C(C)(C)O',
                                      'name': '(R)-(+)-alpha-terpineol',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'},
                                  {   'smiles': 'CC(C)C1(O)CCC(C)=CC1',
                                      'name': '4-terpineol',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'},
                                  {   'smiles': 'Clc1c2COC(=O)c2c(Cl)c(Cl)c1Cl',
                                      'name': '4,5,6,7-tetrachloro-2-benzofuran-1(3H)-one',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'},
                                  {   'smiles': 'CC1=CCC(CC1)C(C)(C)O',
                                      'name': 'alpha-terpineol',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'},
                                  {   'smiles': 'OC1(CCC(CC1)C(C)=C)C',
                                      'name': 'beta-terpineol',
                                      'reason': 'No primary alcohol structure '
                                                'found in SMILES'}],
    'attempt': 1,
    'reasoning': 'Based on the provided information and previous outcomes, it '
                 'seems the SMARTS pattern used to detect primary alcohols was '
                 'insufficient. The pattern may not have captured the full '
                 'range of primary alcohol structures, specifically those that '
                 'don\'t fit the basic `"[CX3H2][OH]"` pattern due to complex '
                 'connectivity or longer alkyl chains which are still attached '
                 'via primary bonding.\n'
                 '\n'
                 'In a primary alcohol, the hydroxyl (-OH) group should be '
                 'directly attached to a carbon that:\n'
                 '1. Is attached to at least two hydrogen atoms (or).\n'
                 '2. Is attached to one other carbon and two hydrogen atoms '
                 '(relax the need for fully saturated cases).\n'
                 '\n'
                 'Here is how we can improve the pattern:\n'
                 '- Expand the SMARTS pattern to cover more complex cases, '
                 'like carbon with many surrounding atoms (e.g. cyclic '
                 'structures, or those bound at terminal carbons with multiple '
                 'attachments).\n'
                 '\n'
                 "Let's revise the code to improve the detection of primary "
                 'alcohols:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 780,
    'num_false_positives': 100,
    'num_true_negatives': 117,
    'num_false_negatives': 9,
    'num_negatives': None,
    'precision': 0.8863636363636364,
    'recall': 0.9885931558935361,
    'f1': 0.9346914319952068,
    'accuracy': 0.8916500994035785,
    'negative_predictive_value': 0.9285714285714286}