"""
Classifies: CHEBI:35179 2-oxo monocarboxylic acid anion
"""
"""
Classifies: CHEBI:38478 2-oxo monocarboxylic acid anion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2_oxo_monocarboxylic_acid_anion(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid anion based on its SMILES string.
    Must have a carboxylate group with an oxo group at the 2-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid anion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylate group pattern ([O-]C(=O)-)
    carboxylate_pattern = Chem.MolFromSmarts("[O-]C(=O)")
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    
    if not carboxylate_matches:
        return False, "No carboxylate group found"
    
    if len(carboxylate_matches) > 1:
        return False, "More than one carboxylate group found"

    # Look for the 2-oxo pattern: [O-]C(=O)C(=O)R
    oxo_acid_pattern = Chem.MolFromSmarts("[O-]C(=O)C(=O)[#6,#7,#8,#9,#16,#15]")
    if not mol.HasSubstructMatch(oxo_acid_pattern):
        return False, "No 2-oxo group adjacent to carboxylate found"
    
    # Make sure there's only one ketone group
    ketone_pattern = Chem.MolFromSmarts("[#6]C(=O)[#6]")
    additional_ketones = mol.GetSubstructMatches(ketone_pattern)
    
    # The 2-oxo group will be counted in additional_ketones, so we check if there are more than 1
    if len(additional_ketones) > 1:
        return False, "Additional ketone groups found beyond the 2-position"
    
    # Additional check to ensure the oxo group is at position 2
    # First get the carboxylate carbon atom
    carboxylate_carbon = carboxylate_matches[0][1]  # Index 1 is the carbon of C(=O)[O-]
    
    # Get all ketone carbons
    ketone_carbons = [match[1] for match in mol.GetSubstructMatches(Chem.MolFromSmarts("C(=O)"))]
    
    # Check if any ketone carbon is directly bonded to the carboxylate carbon
    for ketone_carbon in ketone_carbons:
        bond = mol.GetBondBetweenAtoms(carboxylate_carbon, ketone_carbon)
        if bond is not None:
            return True, "Contains carboxylate with oxo group at 2-position"
    
    return False, "Oxo group not at 2-position"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35179',
                          'name': '2-oxo monocarboxylic acid anion',
                          'definition': 'An oxo monocarboxylic acid anion in '
                                        'which the oxo group is located at the '
                                        '2-position.',
                          'parents': ['CHEBI:35902'],
                          'xrefs': [   'MetaCyc:2-Oxo-carboxylates',
                                       'PMID:10850983'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No 2-oxo group adjacent to '
                                               'carboxylate found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No carboxylate group found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No 2-oxo group adjacent to '
                                               'carboxylate found'}],
    'sample_false_negatives': [   {   'smiles': 'C1=C[C@H](C=C[C@]1(C([O-])=O)CC(C(=O)[O-])=O)[NH3+]',
                                      'name': '4-amino-4-deoxyprephenate',
                                      'reason': 'More than one carboxylate '
                                                'group found'},
                                  {   'smiles': '[NH3+][C@@H](CCCC(=O)C([O-])=O)C([O-])=O',
                                      'name': '(S)-2-amino-6-oxopimelate',
                                      'reason': 'More than one carboxylate '
                                                'group found'},
                                  {   'smiles': 'CC(=O)NCCCCC([O-])=O',
                                      'name': '5-acetamidopentanoate',
                                      'reason': 'No 2-oxo group adjacent to '
                                                'carboxylate found'},
                                  {   'smiles': '[H]C(=O)C([O-])=O',
                                      'name': 'glyoxylate',
                                      'reason': 'No 2-oxo group adjacent to '
                                                'carboxylate found'},
                                  {   'smiles': 'CC(=O)CC(=O)C([O-])=O',
                                      'name': 'acetylpyruvate',
                                      'reason': 'Additional ketone groups '
                                                'found beyond the 2-position'},
                                  {   'smiles': 'CC(=O)C(C)(O)C([O-])=O',
                                      'name': '2-acetyllactate',
                                      'reason': 'No 2-oxo group adjacent to '
                                                'carboxylate found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify 2-oxo '
                 'monocarboxylic acid anions. Looking at the examples, these '
                 'compounds have:\n'
                 '\n'
                 '1. A carboxylate group ([O-]C(=O)-)\n'
                 '2. An oxo (ketone) group at the 2-position (C(=O))\n'
                 '3. These two groups must be adjacent to each other\n'
                 '4. Only one carboxylate group should be present\n'
                 '5. The oxo group must be at position 2 (adjacent to the '
                 'carboxylate)\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 73,
    'num_false_positives': 43,
    'num_true_negatives': 142178,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.6293103448275862,
    'recall': 0.9240506329113924,
    'f1': 0.7487179487179487,
    'accuracy': 0.9996556570625439,
    'negative_predictive_value': 0.9999578011590615}