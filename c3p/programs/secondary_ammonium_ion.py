"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: secondary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule contains a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion has an NH2+ group with exactly two carbons attached.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all nitrogen atoms
    potential_ammonium = False
    reason = "No secondary ammonium ion found"
    
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            # Check formal charge
            if atom.GetFormalCharge() != 1:
                continue
                
            # Count attached hydrogens
            n_hydrogens = atom.GetTotalNumHs()
            
            # Count attached carbons
            carbon_neighbors = sum(1 for neighbor in atom.GetNeighbors() 
                                if neighbor.GetAtomicNum() == 6)
            
            # Secondary ammonium should have:
            # - Exactly 2 hydrogens (NH2+)
            # - Exactly 2 carbon neighbors
            # - Formal charge of +1
            if n_hydrogens == 2 and carbon_neighbors == 2:
                return True, "Contains NH2+ group with exactly two carbons attached"
            
            # Help identify why classification failed
            if n_hydrogens > 2:
                reason = "Primary ammonium ion (NH3+)"
            elif n_hydrogens == 1 and carbon_neighbors == 3:
                reason = "Tertiary ammonium ion (NHR3+)"
            elif n_hydrogens == 0 and carbon_neighbors == 4:
                reason = "Quaternary ammonium ion (NR4+)"
            elif carbon_neighbors < 2:
                reason = f"Only {carbon_neighbors} carbon(s) attached to NH2+ group"
            elif n_hydrogens < 2:
                reason = f"Only {n_hydrogens} hydrogen(s) on charged nitrogen"
                
    return False, reason

__metadata__ = {
    'chemical_class': {
        'name': 'secondary ammonium ion',
        'definition': 'An organic cation obtained by protonation of any secondary amino compound; major species at pH 7.3.',
        'parents': ['organic cation', 'ammonium ion']
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:137419',
                          'name': 'secondary ammonium ion',
                          'definition': 'An organic cation obtained by '
                                        'protonation of any secondary amino '
                                        'compound; major species at pH 7.3.',
                          'parents': ['CHEBI:25697', 'CHEBI:35274'],
                          'xrefs': ['MetaCyc:Secondary-Amines'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O=C1O[C@H]2[C@H](O[C@H](C2)C[C@H](O)C)C=3C1=C(O)C(OC)=C(OC)C3',
                                     'name': '(12R)-12-hydroxymonocerin',
                                     'reason': 'No secondary ammonium ion '
                                               'found'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCCCC/C=C\\CCCCCCCC)(OCCN)(O)=O',
                                     'name': 'PE(20:1(11Z)/16:0)',
                                     'reason': 'No secondary ammonium ion '
                                               'found'},
                                 {   'smiles': 'O=C(C1=C(N)C2=C(OC(C)(C)C=C2)C=C1)C[C@@H](NC(=O)C)CO',
                                     'name': 'Fusarochromene',
                                     'reason': 'No secondary ammonium ion '
                                               'found'},
                                 {   'smiles': 'COc1cc(\\C=C\\C(=O)c2c(O)cc(O)cc2O)ccc1O',
                                     'name': 'homoeriodictyol chalcone',
                                     'reason': 'No secondary ammonium ion '
                                               'found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](O[C@@H](O[C@H]2[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]2CO)O[C@@H]([C@H](O)[C@@H](NC(=O)C)CO)[C@H](O)CO)[C@H]1O)CO[C@H]3O[C@@H]([C@@H](O)[C@H](O)[C@@H]3O[C@@H]4O[C@@H]([C@@H](O[C@@H]5O[C@@H]([C@H](O)[C@H](O[C@]6(O[C@H]([C@H](NC(=O)C)[C@@H](O)C6)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]5O)CO)[C@H](O)[C@H]4NC(=O)C)CO)CO)[C@H]7O[C@@H]([C@@H](O)[C@H](O)[C@@H]7O[C@@H]8O[C@@H]([C@@H](O[C@@H]9O[C@@H]([C@H](O)[C@H](O)[C@H]9NC(=O)C)CO[C@]%10(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%10)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]8NC(=O)C)CO)CO',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[[(2R,3R,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4S,5S,6R)-2-[(2R,3R,4S,5S,6S)-2-[[(2S,3S,4S,5S,6R)-3-[(2S,3R,4R,5S,6R)-3-acetamido-5-[(2S,3R,4S,5S,6R)-4-[(2S,4S,5R,6R)-5-acetamido-2-carboxy-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxan-2-yl]oxy-3,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4-hydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]-6-[(2R,3S,4R,5R,6S)-5-acetamido-6-[(2R,3S,4R,5S)-5-acetamido-1,2,4,6-tetrahydroxyhexan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,5-dihydroxyoxan-4-yl]oxy-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]oxy-4-hydroxy-2-(hydroxymethyl)oxan-3-yl]oxy-3,4-dihydroxyoxan-2-yl]methoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No secondary ammonium ion '
                                               'found'},
                                 {   'smiles': 'O1C=2C(CC=C(C)C)=C(OC)C=C(O)C2C(=O)C=C1C3=CC=C(O)C=C3',
                                     'name': 'Artonin U',
                                     'reason': 'No secondary ammonium ion '
                                               'found'},
                                 {   'smiles': 'O=C(CCCCCCCCC)CCCCCC(O)=O',
                                     'name': '7-Keto palmitic acid',
                                     'reason': 'No secondary ammonium ion '
                                               'found'},
                                 {   'smiles': 'CNC(=N)NCCC[C@H]([NH3+])C([O-])=O',
                                     'name': 'N(omega)-methyl-L-arginine '
                                             'zwitterion',
                                     'reason': 'Primary ammonium ion (NH3+)'},
                                 {   'smiles': '[C@]([C@@](/C=C/CCCCCCCCCCC)(O)[H])(NC(=O)C(CCCCCCCCCCCCCCCCCCCCCC)O)([H])CO',
                                     'name': 'N-(2-hydroxy-tetracosanoyl)-hexadecasphing-4-enine',
                                     'reason': 'No secondary ammonium ion '
                                               'found'},
                                 {   'smiles': 'COc1cc2oc(-c3ccc(O)cc3O)c(CC=C(C)C)c(=O)c2c(O)c1\\C=C\\C(C)C',
                                     'name': 'artocarpin',
                                     'reason': 'No secondary ammonium ion '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'C(CC[NH3+])CC[NH+](C)C',
                                      'name': 'N,N-dimethyl-cadaverine(+2)',
                                      'reason': 'Tertiary ammonium ion '
                                                '(NHR3+)'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify secondary ammonium '
                 "ions. Here's my thought process:\n"
                 '\n'
                 '1. A secondary ammonium ion has:\n'
                 '   - A positively charged nitrogen (NH2+)\n'
                 '   - Two carbon atoms attached to the nitrogen (making it '
                 'secondary)\n'
                 '   - The nitrogen should not have more than 2 carbons '
                 'attached (would be tertiary)\n'
                 '   - Not be a primary ammonium (NH3+) or quaternary ammonium '
                 '(N+)\n'
                 '\n'
                 '2. Key checks needed:\n'
                 '   - Presence of NH2+ group\n'
                 '   - Count of carbon atoms attached to the nitrogen\n'
                 '   - Verify charge state\n'
                 '   - Check for correct valence\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 26,
    'num_false_positives': 100,
    'num_true_negatives': 37071,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.20634920634920634,
    'recall': 0.9629629629629629,
    'f1': 0.3398692810457516,
    'accuracy': 0.9972848002580784,
    'negative_predictive_value': 0.999973025463962}