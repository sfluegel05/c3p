"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
from rdkit import Chem

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic steroid backbone pattern: 17 carbon atoms arranged in four fused rings
    # This example uses a general steroid core scaffold (cyclopenta[a]phenanthrene) as a simple pattern
    steroid_pattern = Chem.MolFromSmarts('C1CCC2C(C1)CCC3C2CCC4C3CCC4') 
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Basic steroid backbone not found"

    # Look for the 11beta-hydroxy group: CCOH at position 11 with beta configuration
    # Simplified SMARTS assuming the general placement
    # Beta configuration commonly signifies a distinct stereochemistry
    hydroxy_11beta_pattern = Chem.MolFromSmarts('[C@@H](O)[C@]1CCCCC1')
    if not mol.HasSubstructMatch(hydroxy_11beta_pattern):
        return False, "No 11beta-hydroxy group found with correct configuration"

    return True, "Contains basic steroid backbone with an 11beta-hydroxy group of the correct configuration"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35346',
                          'name': '11beta-hydroxy steroid',
                          'definition': 'Any 11-hydroxy steroid in which the '
                                        'hydroxy group at position 11 has '
                                        'beta- configuration.',
                          'parents': ['CHEBI:36841'],
                          'xrefs': ['KEGG:C01058'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@@H]1OC(=O)C=C1',
                                     'name': '(S)-5-methylfuran-2(5H)-one',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'Cc1cccc(c1)C(O)=O',
                                     'name': 'm-toluic acid',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'CN(C)CC(=O)NC[C@@H]1[C@@H]([C@@H](N1)CO)C2=CC=C(C=C2)C3=CC=CC(=C3)C#N',
                                     'name': 'N-[[(2S,3S,4R)-3-[4-(3-cyanophenyl)phenyl]-4-(hydroxymethyl)-2-azetidinyl]methyl]-2-(dimethylamino)acetamide',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OCCCCCCCCCC(=O)CC(O)=O)[C@H](O)C[C@H]1O',
                                     'name': 'bkos#20',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'C([C@@](COC(CCCCCCCCC/C=C\\CCCCCC)=O)(OC(CCCCCCCCC/C=C\\CCCCCC)=O)[H])OP(=O)(O)OC[C@@](CO)([H])O',
                                     'name': 'PG(18:1(11Z)/18:1(11Z))',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'O=C(N1C(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(C(=O)NC(CO)CC=2C3=C(C=CC=C3)NC2)CCC(=O)N)(C)C)(C)C)(CC)C)CCC1)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)CNC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C(NC(=O)C)CC4=CC=CC=C4)(C)C)CO)(C)C)(C)C)CC(C)C)CCC(=O)N)(C)C)(C)C)C)C)(C)C',
                                     'name': 'Chrysospermin B',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'CC1(C=CC2=C3C(=CC(=C2O1)OC)C(=O)C(=CO3)C4=CC5=C(C=C4OC)OCO5)C',
                                     'name': '6-methoxy-3-(6-methoxy-1,3-benzodioxol-5-yl)-8,8-dimethyl-4-pyrano[2,3-h][1]benzopyranone',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)NCC(O)=O)C',
                                     'name': 'Glu-Thr-Gly',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'O([C@H]1[C@H](O[C@@H]2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)[C@@H](NC(=O)C)[C@@H](O[C@@H]1CO)O[C@H]3[C@@H](O)[C@H](O)[C@H](O[C@@H]3OC[C@H]4O[C@@H](O[C@H]5[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]5CO)O)[C@@H](O)[C@@H](O)[C@@H]4O)CO)[C@@H]6O[C@@H]([C@H](O)[C@H](O[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]6O)CO',
                                     'name': '(2S,4S,5R,6R)-5-Acetamido-2-[(2S,3R,4S,5S,6R)-2-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2S,3S,4S,5S,6R)-2-[[(2R,3S,4S,5S,6S)-6-[(2R,3S,4R,5R)-5-acetamido-4,6-dihydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4,5-trihydroxyoxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-2-(hydroxymethyl)-4-[(2S,3S,4R,5S,6S)-3,4,5-trihydroxy-6-methyloxan-2-yl]oxyoxan-3-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-4-yl]oxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'Basic steroid backbone not '
                                               'found'},
                                 {   'smiles': 'COCCOC(C)C(O)=O',
                                     'name': '2-(2-methoxyethoxy)propanoic '
                                             'acid',
                                     'reason': 'Basic steroid backbone not '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12C[C@@H](C)[C@H](C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1(Cl)[C@@]2([H])C[C@H](F)C2=CC(=O)C=C[C@]12C',
                                      'name': 'clocortolone',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)COC(=O)CN(CC)CC',
                                      'name': 'hydrocortamate',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': '[H][C@@]12C[C@@H](C)[C@](OC(=O)CC)(C(=O)SCF)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])C[C@H](F)C2=CC(=O)C=C[C@]12C',
                                      'name': 'fluticasone propionate',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@]1(C)C(=O)CC[C@@]21[H]',
                                      'name': '11beta-hydroxyandrost-4-ene-3,17-dione',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': 'CC(=O)OCC(=O)[C@@]1(O)CC[C@H]2[C@@H]3C[C@H](F)C4=CC(=O)CC[C@]4(C)[C@H]3[C@@H](O)C[C@]12C',
                                      'name': '6alpha-fluoro-17-hydroxycorticosterone '
                                              '21-acetate',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(CO)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)CC[C@]12C',
                                      'name': '18-hydroxycortisol',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)CO',
                                      'name': 'cortisol',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)COS(O)(=O)=O',
                                      'name': 'cortisol 21-sulfate',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': '[H][C@@]12C[C@@H](C)[C@](O)(C(=O)SCF)[C@@]1(C)C[C@H](O)[C@@]1(F)[C@@]2([H])C[C@H](F)C2=CC(=O)C=C[C@]12C',
                                      'name': 'fluticasone',
                                      'reason': 'Basic steroid backbone not '
                                                'found'},
                                  {   'smiles': '[H][C@@]12CCc3cc(O)ccc3[C@@]1([H])[C@@H](O)C[C@]1(C)[C@@H](O)CC[C@@]21[H]',
                                      'name': '11beta-hydroxy-17beta-estradiol',
                                      'reason': 'Basic steroid backbone not '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'To classify a chemical entity as an 11beta-hydroxy steroid, '
                 'we need to determine if it contains a hydroxy group '
                 'specifically at position 11 of the steroid backbone, with a '
                 'beta configuration. Steroid backbones typically consist of '
                 'four carbon rings arranged in a specific pattern '
                 '(cyclopenta[a]phenanthrene).\n'
                 '\n'
                 "Here's the step-by-step reasoning for creating a function "
                 '`is_11beta_hydroxy_steroid`:\n'
                 '\n'
                 '1. **Setup RDKit**: Import the necessary modules from RDKit '
                 'for handling chemical structures and substructure '
                 'searching.\n'
                 '2. **Parse SMILES**: Convert the given SMILES string to an '
                 'RDKit molecule object.\n'
                 '3. **Basic Steroid Structure**: Verify the basic structure '
                 'of a steroid, which generally includes four fused carbon '
                 'rings.\n'
                 '4. **11beta-hydroxy Group Identification**: Identify if '
                 "there's a hydroxyl group at the 11th carbon in the steroid "
                 'structure with beta configuration.\n'
                 '5. **Return Result**: Return a boolean indicating if the '
                 'input is an 11beta-hydroxy steroid and a reason for the '
                 'classification.\n'
                 '\n'
                 "Let's put this into code now:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 5,
    'num_false_positives': 100,
    'num_true_negatives': 14746,
    'num_false_negatives': 66,
    'num_negatives': None,
    'precision': 0.047619047619047616,
    'recall': 0.07042253521126761,
    'f1': 0.05681818181818182,
    'accuracy': 0.9888717570557083,
    'negative_predictive_value': 0.9955441533891439}