"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: 6-aminopurines
Definition: Any compound having 6-aminopurine (adenine) as part of its structure
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule contains a 6-aminopurine (adenine) moiety based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains 6-aminopurine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for 6-aminopurine (adenine)
    # This pattern represents the fused bicyclic system with amino group at position 6
    # The pattern matches:
    # - A 5-membered ring fused to a 6-membered ring (purine scaffold)
    # - An amino group (-NH2) at position 6
    # - The correct number and position of nitrogens in the rings
    adenine_pattern = Chem.MolFromSmarts('c1[nH]c2c(n1)c(N)nc[nH]2')
    
    # Alternative pattern that also matches N9-substituted adenines (as in nucleotides)
    adenine_pattern_2 = Chem.MolFromSmarts('c1nc2c(n1)c(N)ncn2')
    
    # Check for matches
    has_adenine_1 = mol.HasSubstructMatch(adenine_pattern)
    has_adenine_2 = mol.HasSubstructMatch(adenine_pattern_2)
    
    if has_adenine_1 or has_adenine_2:
        return True, "Contains 6-aminopurine (adenine) moiety"
    
    return False, "Does not contain 6-aminopurine structure"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:20706',
                          'name': '6-aminopurines',
                          'definition': 'Any compound having 6-aminopurine '
                                        '(adenine) as part of its structure.',
                          'parents': ['CHEBI:22527'],
                          'xrefs': [   'PMID:1646334',
                                       'PMID:18524423',
                                       'PMID:7342604'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'C[C@H]1CN(C(=O)C2=C(C3=CC=CC=C3CO[C@@H]1CN(C)C(=O)NC4=CC=C(C=C4)OC)C5=CC=CC=C5N2C)[C@H](C)CO',
                                     'name': 'LSM-30659',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'c1c[nH+]c[nH]1',
                                     'name': 'imidazolium cation',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'ClC=1C(O)=C2O[C@@]([C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC)[C@@H](O)C(=C2)C1)C(C)=C)C(=O)N3[C@@H](C=CC3)C(=O)N/C(=C(/CC)\\C)/C(=O)N/C(=C/C(O)=O)/C(O)=O)(CC)C',
                                     'name': 'Phomopsin A',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'ClC(Cl)[C@H](O)CC=1OC(=O)C=2C(O)=CC(=CC2C1)O',
                                     'name': 'Desmethyldichlorodiaportin',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'ClC=1C(=C(O)C2=C(C1O)C(=O)C=CC2=O)CC=C(C)C',
                                     'name': 'Chlorosesamone',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'CCCCC[C@H]1O[C@H]1C\\C=C/CCCCCCCC(O)=O',
                                     'name': '(+)-vernolic acid',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'O(C1C(O)C(OC(OC=2C3=C(C=CC2)C=C(C(=C3O)C(=O)C)C)C1O)CO)C4OC(C(O)C(O)C4O)CO',
                                     'name': 'Orientaloside',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'O(C1=C(C=CC=C1C(O)=O)C)C(=O)C',
                                     'name': 'CRESOPYRINE',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'C[C@H]1CN(S(=O)(=O)C2=C(C=C(C=C2)C#CC3=CC=NC=C3)O[C@@H]1CN(C)C(=O)NC4=CC5=C(C=C4)OCO5)[C@@H](C)CO',
                                     'name': '3-(1,3-benzodioxol-5-yl)-1-[[(4S,5S)-2-[(2S)-1-hydroxypropan-2-yl]-4-methyl-1,1-dioxo-8-(2-pyridin-4-ylethynyl)-4,5-dihydro-3H-6,1$l^{6},2-benzoxathiazocin-5-yl]methyl]-1-methylurea',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'},
                                 {   'smiles': 'O1C2C3C(CCC3=C)C(CCC2C(C1=O)=C)=C',
                                     'name': '3,6,9-Trimethylidene-3a,4,5,6a,7,8,9a,9b-octahydroazuleno[4,5-b]furan-2-one',
                                     'reason': 'Does not contain 6-aminopurine '
                                               'structure'}],
    'sample_false_negatives': [   {   'smiles': 'O=NN1CCNCC1',
                                      'name': 'N-Mononitrosopiperazine',
                                      'reason': 'Does not contain '
                                                '6-aminopurine structure'},
                                  {   'smiles': 'O1C(CCC1)C(=O)N2CCNCC2',
                                      'name': '1-(Tetrahydro-2-furoyl)piperazine',
                                      'reason': 'Does not contain '
                                                '6-aminopurine structure'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify 6-aminopurines. '
                 'Looking at the examples, I notice that all these compounds '
                 'contain the adenine (6-aminopurine) structure, which '
                 'consists of a purine ring system with an amino group at '
                 'position 6.\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 676,
    'num_false_positives': 100,
    'num_true_negatives': 9180,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.8711340206185567,
    'recall': 0.9970501474926253,
    'f1': 0.9298486932599724,
    'accuracy': 0.9897569793131151,
    'negative_predictive_value': 0.999782182531039}