"""
Classifies: CHEBI:22307 aldoxime
"""
"""
Classifies: CHEBI:33566 aldoxime
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aldoxime(smiles: str):
    """
    Determines if a molecule is an aldoxime based on its SMILES string.
    Aldoximes are compounds with the general structure R-CH=N-OH, where R is any carbon group.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an aldoxime, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for oxime group pattern (C=N-OH)
    # [CX3H1] ensures carbon has exactly one hydrogen
    # [NX2] ensures nitrogen has double bond
    # [OX2H1] ensures oxygen is connected to one hydrogen
    oxime_pattern = Chem.MolFromSmarts("[CX3H1]=[NX2]-[OX2H1]")
    
    # Alternative pattern for deprotonated oximes
    oxime_anion_pattern = Chem.MolFromSmarts("[CX3H1]=[NX2]-[OX1-]")
    
    matches = mol.GetSubstructMatches(oxime_pattern)
    anion_matches = mol.GetSubstructMatches(oxime_anion_pattern)
    
    total_matches = len(matches) + len(anion_matches)
    
    if total_matches == 0:
        return False, "No aldoxime group (C(H)=N-OH) found"
    
    # Check if the carbon of C=N is connected to exactly one hydrogen
    # This distinguishes aldoximes from ketoximes
    for match in matches:
        carbon_atom = mol.GetAtomWithIdx(match[0])
        if carbon_atom.GetTotalNumHs() != 1:
            return False, "Oxime carbon has incorrect number of hydrogens"
            
    for match in anion_matches:
        carbon_atom = mol.GetAtomWithIdx(match[0])
        if carbon_atom.GetTotalNumHs() != 1:
            return False, "Oxime carbon has incorrect number of hydrogens"
    
    # Verify carbon is not part of C=O group (would indicate an oxime ester)
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=O)-[OX2]")
    if mol.HasSubstructMatch(carbonyl_pattern):
        # Check if the carbonyl carbon is the same as oxime carbon
        carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
        for oxime_match in matches:
            oxime_carbon = oxime_match[0]
            for carbonyl_match in carbonyl_matches:
                if oxime_carbon == carbonyl_match[0]:
                    return False, "Contains oxime ester group instead of aldoxime"
    
    return True, "Contains aldoxime group (R-CH=N-OH)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:22307',
                          'name': 'aldoxime',
                          'definition': 'Oximes of aldehydes RCH=NOH.',
                          'parents': ['CHEBI:25750'],
                          'xrefs': ['KEGG:C02658'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)NC4=CC=C(C=C4)OC)C5=CC=CC=C5N2C)[C@H](C)CO',
                                     'name': 'LSM-30659',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'ClC=1C(O)=C2O[C@@]([C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC)[C@@H](O)C(=C2)C1)C(C)=C)C(=O)N3[C@@H](C=CC3)C(=O)N/C(=C(/CC)\\C)/C(=O)N/C(=C/C(O)=O)/C(O)=O)(CC)C',
                                     'name': 'Phomopsin A',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'CCCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O',
                                     'name': '(+)-vernolic acid',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'O(C1C(O)C(OC(OC=2C3=C(C=CC2)C=C(C(=C3O)C(=O)C)C)C1O)CO)C4OC(C(O)C(O)C4O)CO',
                                     'name': 'Orientaloside',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=NC=C3)O[C@@H]1CN(C)C(=O)NC4=CC5=C(C=C4)OCO5)[C@@H](C)CO',
                                     'name': '3-(1,3-benzodioxol-5-yl)-1-[[(4S,5S)-2-[(2S)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-8-(2-pyridin-4-ylethynyl)-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-1-methylurea',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'},
                                 {   'smiles': 'O1C2C3C(CCC3=C)C(CCC2C(C1=O)=C)=C',
                                     'name': '3,6,9-Trimethylidene-3a,4,5,6a,7,8,9a,9b-octahydroazuleno[4,5-b]furan-2-one',
                                     'reason': 'No aldoxime group (C(H)=N-OH) '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'CSC(=NO)C(=O)N(C)C',
                                      'name': '2-(hydroxyimino)-N,N-dimethyl-2-(methylsulfanyl)acetamide',
                                      'reason': 'No aldoxime group (C(H)=N-OH) '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'Let me create a program to classify aldoximes. The key '
                 'characteristics of aldoximes are:\n'
                 '\n'
                 '1. They contain a C=N-OH group (oxime group)\n'
                 '2. The carbon of the C=N group has one hydrogen attached '
                 '(aldehyde-derived)\n'
                 '3. The other substituent on the C=N can be any carbon group '
                 '(R)\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 48,
    'num_false_positives': 28,
    'num_true_negatives': 142223,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.631578947368421,
    'recall': 0.9795918367346939,
    'f1': 0.768,
    'accuracy': 0.9997962052002811,
    'negative_predictive_value': 0.9999929688378896}