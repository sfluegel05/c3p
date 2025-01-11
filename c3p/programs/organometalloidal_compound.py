"""
Classifies: CHEBI:143084 organometalloidal compound
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_organometalloidal_compound(smiles: str):
    """
    Determines if a molecule is an organometalloidal compound based on its SMILES string.
    An organometalloidal compound has bonds between metalloid atoms (primarily As, but also B, Si, Ge, Sb, Te)
    and carbon atoms of an organyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organometalloidal compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define metalloids (focusing primarily on As since it's most common in examples)
    metalloids = [33]  # As atomic number
    
    # Track if we found valid metalloid-carbon bonds
    found_metalloid_carbon = False
    
    # Examine each atom in the molecule
    for atom in mol.GetAtoms():
        # Check if atom is a metalloid
        if atom.GetAtomicNum() in metalloids:
            # Look for carbon neighbors
            carbon_neighbors = [n for n in atom.GetNeighbors() 
                              if n.GetAtomicNum() == 6]
            
            if carbon_neighbors:
                # Check if any carbon is part of an organyl group
                for carbon in carbon_neighbors:
                    # Skip if carbon is part of CO3 or similar inorganic groups
                    oxygen_count = sum(1 for n in carbon.GetNeighbors() 
                                     if n.GetAtomicNum() == 8)
                    if oxygen_count >= 3:
                        continue
                        
                    # Valid metalloid-carbon bond found
                    found_metalloid_carbon = True
                    break
    
    if not found_metalloid_carbon:
        return False, "No valid metalloid-carbon bonds found"
        
    # Additional checks for arsenic compounds
    as_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 33]
    if not as_atoms:
        return False, "No arsenic atoms found"
        
    # Count total number of As-C bonds
    total_as_c_bonds = sum(len([n for n in as_atom.GetNeighbors() 
                               if n.GetAtomicNum() == 6]) 
                          for as_atom in as_atoms)
    
    # Verify that we have at least one As-C bond
    if total_as_c_bonds == 0:
        return False, "No As-C bonds found"
        
    # Success case
    return True, f"Found {total_as_c_bonds} As-C bond(s) in organometalloidal compound"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:143084',
                          'name': 'organometalloidal compound',
                          'definition': 'A compound having bonds between one '
                                        'or more metalloid atoms and one or '
                                        'more carbon atoms of an organyl '
                                        'group.',
                          'parents': ['CHEBI:50860'],
                          'xrefs': [   'PMID:15246002',
                                       'PMID:17741021',
                                       'PMID:8403081',
                                       'PMID:908316'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: unterminated string literal (detected at line 9) '
               '(<string>, line 9)\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@H]2[C@H](O[C@H](C2)C[C@H](O)C)C=3C1=C(O)C(OC)=C(OC)C3',
                                     'name': '(12R)-12-hydroxymonocerin',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                     'name': 'PE(20:1(11Z)/16:0)',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': 'O=C(C1=C(N)C2=C(OC(C)(C)C=C2)C=C1)C[C@@H](NC(=O)C)CO',
                                     'name': 'Fusarochromene',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': 'COc1cc(\\C=C\\C(=O)c2c(O)cc(O)cc2O)ccc1O',
                                     'name': 'homoeriodictyol chalcone',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9NC(=O)C)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5S)-5-acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': 'O1C=2C(CC=C(C)C)=C(OC)C=C(O)C2C(=O)C=C1C3=CC=C(O)C=C3',
                                     'name': 'Artonin U',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': 'O=C(CCCCCCCCC)CCCCCC(O)=O',
                                     'name': '7-Keto palmitic acid',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': 'CNC(=N)NCCC[C@H]([NH3+])C([O-])=O',
                                     'name': 'N(omega)-methyl-L-arginine '
                                             'zwitterion',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': '[C@]([C@@](/C=C/CCCCCCCCCCC)(O)[H])(NC(=O)C(CCCCCCCCCCCCCCCCCCCCCC)O)([H])CO',
                                     'name': 'N-(2-hydroxy-tetracosanoyl)-hexadecasphing-4-enine',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'},
                                 {   'smiles': 'COc1cc2oc(-c3ccc(O)cc3O)c(CC=C(C)C)c(=O)c2c(O)c1\\C=C\\C(C)C',
                                     'name': 'artocarpin',
                                     'reason': 'No valid metalloid-carbon '
                                               'bonds found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'Let me help analyze the issues and create a new solution.\n'
                 '\n'
                 'Analyzing the examples provided:\n'
                 '1. Most examples contain arsenic (As) as the metalloid\n'
                 '2. The As-C bonds appear in various contexts:\n'
                 '   - Direct methyl groups (C[As])\n'
                 '   - Aromatic rings (c-[As])\n'
                 '   - Aliphatic chains (CC[As])\n'
                 '3. The arsenic can be in different oxidation states and '
                 'forms:\n'
                 '   - Arsonic acids [As](O)(O)=O\n'
                 '   - Arsine groups [As](C)(C)\n'
                 '   - Arsinoyl groups [As](=O)(C)(C)\n'
                 '\n'
                 "Here's a new implementation that should better handle these "
                 'cases:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 29,
    'num_false_positives': 49,
    'num_true_negatives': 142222,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.3717948717948718,
    'recall': 1.0,
    'f1': 0.5420560747663552,
    'accuracy': 0.9996556570625439,
    'negative_predictive_value': 1.0}