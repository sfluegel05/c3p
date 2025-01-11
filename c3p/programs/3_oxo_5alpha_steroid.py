"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
"""
Classifies: 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # First check for basic steroid core (four fused rings)
    steroid_core = Chem.MolFromSmarts("[C]1[C][C][C]2[C]1[C][C][C]3[C]2[C][C][C]4[C][C][C][C][C]34")
    if not mol.HasSubstructMatch(steroid_core):
        return False, "No steroid core structure found"

    # Check for ketone at position 3
    # The pattern looks for a carbonyl (C=O) group in the A ring of the steroid
    oxo_pattern = Chem.MolFromSmarts("[CH2][CH2][C](=[O])[CH2][C](@[H])")
    if not mol.HasSubstructMatch(oxo_pattern):
        return False, "No ketone group at position 3"

    # Check for 5-alpha configuration
    # In 5-alpha steroids, the hydrogen at position 5 is below the plane (alpha)
    # This means the connection between rings A/B is trans
    # Looking for the specific connectivity pattern around C5
    alpha_5_pattern = Chem.MolFromSmarts("[C]1[CH2][C](=O)[CH2][C@@H]([C]2)[CH2][CH2]")
    if not mol.HasSubstructMatch(alpha_5_pattern):
        return False, "No 5-alpha configuration found"

    # Additional check for reasonable molecular weight range of steroids
    mol_wt = sum([atom.GetMass() for atom in mol.GetAtoms()])
    if mol_wt < 250 or mol_wt > 1000:
        return False, "Molecular weight outside typical range for steroids"

    # Count carbons (steroids typically have 19+ carbons)
    carbon_count = len([atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6])
    if carbon_count < 19:
        return False, "Too few carbons for a steroid structure"

    return True, "Molecule contains 3-oxo-5alpha-steroid structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:13601',
                          'name': '3-oxo-5alpha-steroid',
                          'definition': 'A 3-oxo steroid that has alpha '
                                        'configuration at position 5.',
                          'parents': ['CHEBI:47788'],
                          'xrefs': [   'KEGG:C02940',
                                       'MetaCyc:3-Oxo-5-Alpha-Steroids'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@H]2[C@H](O[C@H](C2)C[C@H](O)C)C=3C1=C(O)C(OC)=C(OC)C3',
                                     'name': '(12R)-12-hydroxymonocerin',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                     'name': 'PE(20:1(11Z)/16:0)',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O=C(C1=C(N)C2=C(OC(C)(C)C=C2)C=C1)C[C@@H](NC(=O)C)CO',
                                     'name': 'Fusarochromene',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'COc1cc(\\C=C\\C(=O)c2c(O)cc(O)cc2O)ccc1O',
                                     'name': 'homoeriodictyol chalcone',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9NC(=O)C)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5S)-5-acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O1C=2C(CC=C(C)C)=C(OC)C=C(O)C2C(=O)C=C1C3=CC=C(O)C=C3',
                                     'name': 'Artonin U',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'O=C(CCCCCCCCC)CCCCCC(O)=O',
                                     'name': '7-Keto palmitic acid',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'CNC(=N)NCCC[C@H]([NH3+])C([O-])=O',
                                     'name': 'N(omega)-methyl-L-arginine '
                                             'zwitterion',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': '[C@]([C@@](/C=C/CCCCCCCCCCC)(O)[H])(NC(=O)C(CCCCCCCCCCCCCCCCCCCCCC)O)([H])CO',
                                     'name': 'N-(2-hydroxy-tetracosanoyl)-hexadecasphing-4-enine',
                                     'reason': 'No steroid core structure '
                                               'found'},
                                 {   'smiles': 'COc1cc2oc(-c3ccc(O)cc3O)c(CC=C(C)C)c(=O)c2c(O)c1\\C=C\\C(C)C',
                                     'name': 'artocarpin',
                                     'reason': 'No steroid core structure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'C[C@H](CCCC(C)C(O)=O)[C@H]1CC[C@H]2[C@@H]3CC[C@H]4CC(=O)CC[C@]4(C)[C@H]3CC[C@]12C',
                                      'name': '(5alpha)-3-oxocholestan-26-oic '
                                              'acid',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': 'C[C@@H]1C[C@@]2(C)[C@@H](CC[C@H]3[C@@H]4CC[C@H](OC5CCCCO5)[C@@]4(C)CC[C@H]23)CC1=O',
                                      'name': '2alpha-Methyl-17beta-[(tetrahydro-2H-pyran-2-yl)oxy]-5alpha-androstan-3-one',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': '[C@@]12([C@@]([C@@]3([C@]([C@@H](C(CC3)=O)C)(CC1)[H])C)(CC[C@@]4([C@H](CC[C@@]24[H])O)C)[H])[H]',
                                      'name': '17beta-hydroxy-4alpha-methyl-5alpha-androstan-3-one',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': 'C[C@]1(O)CC[C@H]2[C@@H]3CC[C@H]4CC(=O)CC[C@]4(C)[C@H]3C(=O)C[C@]12C',
                                      'name': '17beta-Hydroxy-17-methyl-5alpha-androstane-3,11-dione',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': 'COC(C)(C)\\C=C\\C[C@@H](C)[C@H]1CC[C@@]2(C)[C@@H]3[C@@H](O)C[C@@H]4[C@]5(C[C@@]35CC[C@]12C)CCC(=O)C4(C)C',
                                      'name': 'combretanone G',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': 'C1[C@]2([C@]3([C@@]([C@@]4([C@](CC([C@@H](C4)C)=O)(CC3)[H])C)(CC[C@@]2([C@H](C1)O)C)[H])[H])[H]',
                                      'name': 'metholone',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': 'C[C@]1(O)CC[C@H]2[C@@H]3CC[C@H]4CC(=O)CC[C@]4(C)[C@@]3(F)[C@@H](O)C[C@]12C',
                                      'name': '11beta,17beta-Dihydroxy-9alpha-fluoro-17alpha-methyl-5alpha-androstan-3-one',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': 'C[C@]12CC[C@H]3[C@@H](CC[C@H]4CC(=O)CC[C@]34C)[C@@H]1CC[C@@H]2OC(=O)c1ccccc1',
                                      'name': 'Stanolone benzoate',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': '[H][C@@]12CC[C@@]3([H])[C@]4([H])CC[C@](O)(C(=O)CO)[C@@]4(C)CC(=O)[C@]3([H])[C@@]1(C)CCC(=O)C2',
                                      'name': '4,5alpha-dihydrocortisone',
                                      'reason': 'No ketone group at position '
                                                '3'},
                                  {   'smiles': 'C[C@@]12[C@@]3([C@]([C@]4([C@@](CC3)(C)[C@@](CC4)(C)O)[H])(CC[C@]1(CC(CC2)=O)[H])[H])[H]',
                                      'name': 'Mestanolone',
                                      'reason': 'No ketone group at position '
                                                '3'}],
    'attempt': 0,
    'reasoning': "I'll create a program to classify 3-oxo-5alpha-steroids. "
                 'These compounds are characterized by:\n'
                 '1. A steroid core structure\n'
                 '2. An alpha (below plane) configuration at position 5\n'
                 '3. A ketone group (=O) at position 3\n'
                 '\n'
                 'Let me write a program that checks for these structural '
                 'features:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 0,
    'num_false_positives': 0,
    'num_true_negatives': 142248,
    'num_false_negatives': 52,
    'num_negatives': None,
    'precision': 0.0,
    'recall': 0.0,
    'f1': 0.0,
    'accuracy': 0.9996345748418833,
    'negative_predictive_value': 0.9996345748418833}