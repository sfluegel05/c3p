"""
Classifies: CHEBI:35627 beta-lactam
"""
from rdkit import Chem

def is_beta_lactam(smiles: str):
    """
    Determines if a molecule is a beta-lactam based on its SMILES string.
    A beta-lactam contains a four-membered lactam (cyclic amide) ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-lactam, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for a four-membered beta-lactam ring
    beta_lactam_pattern = Chem.MolFromSmarts("C1C(=O)NC1")  # Four-membered ring with C=O and N
    
    if mol.HasSubstructMatch(beta_lactam_pattern):
        return True, "Contains a four-membered beta-lactam ring"
    
    return False, "Does not contain a four-membered beta-lactam ring"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35627',
                          'name': 'beta-lactam',
                          'definition': 'A lactam in which the amide bond is '
                                        'contained within a four-membered '
                                        'ring, which includes the amide '
                                        'nitrogen and the carbonyl carbon.',
                          'parents': ['CHEBI:24995'],
                          'xrefs': ['KEGG:C01866', 'Wikipedia:Beta-lactam'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'Does not contain a '
                                               'four-membered beta-lactam '
                                               'ring'}],
    'sample_false_negatives': [   {   'smiles': '[H]C(N[C@@H](C([O-])=O)C(C)(C)S)=C1N=C(Cc2ccccc2)OC1=O',
                                      'name': 'benzylpenicillenate',
                                      'reason': 'Does not contain a '
                                                'four-membered beta-lactam '
                                                'ring'},
                                  {   'smiles': '[H][C@]1(N[C@@H](C([O-])=O)C(C)(C)S1)[C@H](NC(=O)CC1=CC=CC=C1)C(=O)NCCCC',
                                      'name': 'benzylpenicilloyl-butylamine(1-)',
                                      'reason': 'Does not contain a '
                                                'four-membered beta-lactam '
                                                'ring'}],
    'attempt': 0,
    'reasoning': 'Determining whether a compound is a beta-lactam based on its '
                 'SMILES string involves identifying the presence of a '
                 'characteristic four-membered lactam ring that includes the '
                 'amide nitrogen and the carbonyl carbon. The main steps in '
                 'this classification are:\n'
                 '\n'
                 '1. **Parsing the SMILES string**: Utilize RDKit to parse and '
                 'convert the SMILES string into a molecular representation.\n'
                 '2. **SMARTS pattern matching**: Look for the beta-lactam '
                 'substructure, specifically the four-membered ring structure '
                 'with the carbonyl carbon and the nitrogen.\n'
                 '3. **Return classification**: Based on the presence of the '
                 'required structural features, determine if the compound is a '
                 'beta-lactam.\n'
                 '\n'
                 "Let's construct the program following these steps.",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 264,
    'num_false_positives': 100,
    'num_true_negatives': 91158,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.7252747252747253,
    'recall': 0.9924812030075187,
    'f1': 0.8380952380952381,
    'accuracy': 0.9988855382194834,
    'negative_predictive_value': 0.9999780605528741}