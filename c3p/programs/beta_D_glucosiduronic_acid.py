"""
Classifies: CHEBI:15341 beta-D-glucosiduronic acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_beta_D_glucosiduronic_acid(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronic acid based on its SMILES string.
    A beta-D-glucosiduronic acid is a molecule resulting from the formal condensation of any substance with beta-D-glucuronic acid to form a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the beta-D-glucuronic acid moiety SMARTS pattern
    # The anomeric carbon (C1) connected via O-glycosidic bond to any group ([!H])
    # Stereochemistry is specified using @@ and @ symbols
    # The carboxylic acid group at C6 can be protonated or deprotonated ([O-,OH])
    glucuronic_acid_smarts = '[C@@H]1([O][!H])[C@H](O)[C@@H](O)[C@H](O)[C@@H](O1)C(=O)[O-,OH]'
    
    # Create the glucuronic acid pattern molecule
    glucuronic_acid_pattern = Chem.MolFromSmarts(glucuronic_acid_smarts)
    if glucuronic_acid_pattern is None:
        return False, "Failed to create beta-D-glucuronic acid pattern"

    # Check for the beta-D-glucuronic acid substructure
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    return True, "Contains beta-D-glucuronic acid moiety linked via glycosidic bond"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15341',
                          'name': 'beta-D-glucosiduronic acid',
                          'definition': 'A glucosiduronic acid resulting from '
                                        'the formal condensation of any '
                                        'substance with beta-D-glucuronic acid '
                                        'to form a glycosidic bond.',
                          'parents': ['CHEBI:24302'],
                          'xrefs': ['KEGG:C03033'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O[C@H]1[C@@H]([C@@H](C(=O)C1)/C=C/[C@H](O)C/C=C\\C/C=C/C/C=C\\CC)CCC(O)=O',
                                     'name': 'ent-11-D4t-NeuroP',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'C1[C@@H]2[C@@H]([C@@H](N2CC3=CC=CC=C3Cl)CO)C4=CC=CC=C4N1S(=O)(=O)C5=CC=C(C=C5)F',
                                     'name': '[(1R,2aS,8bS)-2-[(2-chlorophenyl)methyl]-4-(4-fluorophenyl)sulfonyl-1,2a,3,8b-tetrahydroazeto[2,3-c]quinolin-1-yl]methanol',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'C=CC1C2C[C@]3([C@@]45C(C6=CC=CC=C6N4)[C@@](CN3)(O[C@@](C25C(=O)OC)(O[C@]1(O[C@@H]7[C@H]([C@@H]([C@H]([C@H](CO)O7)O)O)O)[H])[H])[H])[H]',
                                     'name': 'Cymoside',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'CC1(CC2=C(C(NC3=C2C4=CC=CC=C4C=C3)C5=CC(=C(C=C5)OCC(=O)O)OC)C(=O)C1)C',
                                     'name': '2-[4-(2,2-dimethyl-4-oxo-1,3,5,6-tetrahydrobenzo[a]phenanthridin-5-yl)-2-methoxyphenoxy]acetic '
                                             'acid',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'CC1(C)C2CC(=O)C1(C)CC2O',
                                     'name': '5-hydroxycamphor',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'O1C(C(=O)N2C(CCC2)C(=O)NC(C(CC)C)C(=O)N(C(C(C)C)C(=O)N(C(C(=O)NCCC1=O)C)C)C)CC(O)(C)C',
                                     'name': 'Hydroxydestruxin B',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'CCC1=NN(C(=O)C2=CC3=C(N21)C=CO3)C(C)C(=O)NC(C)CCC4=CC=CO4',
                                     'name': '2-(1-ethyl-4-oxo-3-furo[3,4]pyrrolo[3,5-c][1,2,4]triazinyl)-N-[4-(2-furanyl)butan-2-yl]propanamide',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'O[C@@H]1C(=C)[C@]2(C(C)(C)[C@H](C1)O)CC=C(CO)CC2',
                                     'name': 'Acaciicolinol D',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'C[N+]1(CCCC(C1)NC(=O)CCCC(O)=O)Cc1ccc(cc1)C(=O)NCCO',
                                     'name': '3-(4-carboxybutanamido)-1-{4-[(2-hydroxyethyl)carbamoyl]benzyl}-1-methylpiperidinium',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'},
                                 {   'smiles': 'COC1=CC=CC=C1CNC(=O)C[C@@H]2C[C@@H]3[C@H]([C@@H](O2)CO)OC4=C3C=C(C=C4)NS(=O)(=O)C',
                                     'name': '2-[(1S,3S,4aS,9aR)-1-(hydroxymethyl)-6-(methanesulfonamido)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-[(2-methoxyphenyl)methyl]acetamide',
                                     'reason': 'No beta-D-glucuronic acid '
                                               'moiety found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@]1(C)[C@H](CC[C@@]21[H])O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O',
                                      'name': 'testosterone 17-glucosiduronic '
                                              'acid',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': 'C1[C@]2([C@]3([C@@](C4=C(C=C(O)C=C4)CC3)(CC[C@@]2([C@@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C(O)=O)C)[H])[H])[H]',
                                      'name': 'estriol '
                                              '16-O-(beta-D-glucuronide)',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': 'CC1(C)CC[C@@]2([C@H](O)C[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O[C@@H]7OC[C@H](O)[C@H](O)[C@H]7O[C@@H]7OC[C@@H](O)[C@H](O)[C@H]7O)[C@H]6O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)C(O)=O)C(C)(C)[C@@H]5CC[C@@]34C)[C@@H]2C1)C(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O',
                                      'name': 'gordonoside N',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': 'C1[C@]2([C@]3([C@@](C4=C(C=C(O)C=C4)CC3)(CC[C@@]2([C@H](O[C@@H]5O[C@@H]([C@H]([C@@H]([C@H]5O)O)O)C(O)=O)C1)C)[H])[H])[H]',
                                      'name': '17alpha-estradiol '
                                              '17-O-(beta-D-glucuronide)',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': '[H][C@@]12CC(C)(C)CC(=O)[C@]1(C)CC[C@]1(C)C2=CC[C@]2([H])[C@@]3(C)CC[C@H](O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)[C@](C)(CO)[C@]3([H])CC[C@@]12C',
                                      'name': 'soyasapogenol E '
                                              '3-O-beta-glucuronide',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': 'CC\\C(c1ccccc1)=C(/c1ccccc1)c1ccc(OCC[N+](C)(C)[C@@H]2O[C@@H]([C@@H](O)[C@H](O)[C@H]2O)C(O)=O)cc1',
                                      'name': 'tamoxifen '
                                              'N-beta-D-glucosiduronic acid',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': '[H][C@@]12CC[C@H](O)[C@@]1(C)CC[C@@]1([H])[C@@]2([H])CC[C@@]2([H])C[C@@H](CC[C@]12C)O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(O)=O',
                                      'name': '5alpha-androstane-3beta,17beta-diol '
                                              '3-O-(beta-D-glucuronide)',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': 'C\\C=C(\\C)C(=O)O[C@H]1CC(C)(C)C[C@H]2C3=CC[C@@H]4[C@@]5(C)CC[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O[C@@H]7OC[C@H](O)[C@H](O)[C@H]7O[C@@H]7OC[C@@H](O)[C@H](O)[C@H]7O)[C@H]6O[C@@H]6O[C@H](CO)[C@H](O)[C@H](O)[C@H]6O)C(O)=O)C(C)(C)[C@@H]5CC[C@@]4(C)[C@]3(C)[C@@H](O)[C@@H](O)[C@@]12CO',
                                      'name': 'gordonoside J',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': '[H][C@@]12CC(C)(C)[C@@H](O)[C@@H](O)[C@]1(C)CC[C@]1(C)C2=CC[C@]2([H])[C@@]3(C)CC[C@H](O[C@@H]4O[C@@H]([C@@H](O)[C@H](O)[C@H]4O)C(O)=O)[C@](C)(CO)[C@]3([H])CC[C@@]12C',
                                      'name': 'soyasapogenol A '
                                              '3-O-beta-glucuronide',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'},
                                  {   'smiles': 'C=1(C(CCC(C1C)O[C@H]2[C@@H]([C@H]([C@@H]([C@H](O2)C(O)=O)O)O)O)(C)C)/C=C/C(=C/C=C/C(=C/C(O)=O)/C)/C',
                                      'name': '4-(beta-D-glucopyranuronosyloxy)retinoic '
                                              'acid',
                                      'reason': 'No beta-D-glucuronic acid '
                                                'moiety found'}],
    'attempt': 0,
    'reasoning': 'To determine if a molecule is a beta-D-glucosiduronic acid, '
                 "we'll check for the presence of the beta-D-glucuronic acid "
                 'moiety linked via a glycosidic bond to any other substance. '
                 'The beta-D-glucuronic acid moiety has specific '
                 "stereochemistry and a carboxylic acid group at C6. We'll use "
                 'RDKit to parse the SMILES string and search for a '
                 'substructure pattern that represents this moiety, ensuring '
                 'the anomeric carbon (C1) is connected via an O-glycosidic '
                 'bond.\n'
                 '\n'
                 "Here's the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 83,
    'num_false_positives': 100,
    'num_true_negatives': 22159,
    'num_false_negatives': 32,
    'num_negatives': None,
    'precision': 0.453551912568306,
    'recall': 0.7217391304347827,
    'f1': 0.5570469798657718,
    'accuracy': 0.9941002949852508,
    'negative_predictive_value': 0.9985579739534045}