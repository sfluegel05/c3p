"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: decanoate ester
Definition: A fatty acid ester resulting from the formal condensation of the carboxy group 
of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for ester group pattern (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"
    
    # Pattern for decanoyl group (10-carbon chain with carbonyl)
    # Note: The pattern looks for C(=O)-CCCCCCCCC
    decanoyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4]")
    decanoyl_matches = mol.GetSubstructMatches(decanoyl_pattern)
    
    if not decanoyl_matches:
        return False, "No decanoyl group found"
    
    # Check if any decanoyl group is part of an ester
    found_decanoate_ester = False
    for ester_match in ester_matches:
        for decanoyl_match in decanoyl_matches:
            # Check if the carbonyl carbon of the ester matches the carbonyl carbon of the decanoyl
            if ester_match[1] == decanoyl_match[0]:
                found_decanoate_ester = True
                break
        if found_decanoate_ester:
            break
            
    if not found_decanoate_ester:
        return False, "Decanoyl group not connected via ester linkage"
    
    return True, "Contains decanoyl group connected via ester linkage"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:87658',
                          'name': 'decanoate ester',
                          'definition': 'A fatty acid ester resulting from the '
                                        'formal condensation of the carboxy '
                                        'group of decanoic acid (capric acid) '
                                        'with the hydroxy group of an alcohol '
                                        'or phenol.',
                          'parents': ['CHEBI:35748'],
                          'xrefs': ['PMID:23383323'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'O1[C@H](CC2=C(C=3C=4C(C(O)=C5C3C[C@H](C)OC5)=C(OC)C=C(C4)OC)C6=CC(OC)=CC(=C6C(=C2C1)O)OC)C',
                                     'name': '(3S)-5-[(3S)-10-hydroxy-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-5-yl]-7,9-dimethoxy-3-methyl-3,4-dihydro-1H-benzo[g]isochromen-10-ol',
                                     'reason': 'No ester group found'},
                                 {   'smiles': 'C[C@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)C)C(=O)N(C[C@@H]1OC)C)C)C(=O)C3CCC3',
                                     'name': 'N-[(5R,6S,9R)-8-[cyclobutyl(oxo)methyl]-5-methoxy-3,6,9-trimethyl-2-oxo-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No ester group found'},
                                 {   'smiles': '[C@@]12(CCCC[C@@]1(C=C[C@@H]([C@@H]2CC[C@H](C[C@H](CC(=O)[O-])O)O)C)[H])[H]',
                                     'name': '4a,5-dihydro-ML-236C carboxylate',
                                     'reason': 'No ester group found'},
                                 {   'smiles': 'S(=O)(CC1OC(N2C3=NC=NC(N)=C3N=C2)C(O)C1O)C',
                                     'name': "(S)-5'-Deoxy-5'-(methylsulfinyl)adenosine",
                                     'reason': 'No ester group found'},
                                 {   'smiles': 'OC=1C(=C(C=2C=3C(NC2)=CC=CC3)C(=O)C(=O)C1C=4C=5C(NC4)=CC=CC5)CC=C(C)C',
                                     'name': 'Ochrindole D',
                                     'reason': 'No ester group found'},
                                 {   'smiles': 'O(C1C(O)C(OC1OC2=C(OC=3C(C2=O)=C(O)C=C(OC4OC(C(O)C(O)C4O)C)C3)C5=CC=C(O)C=C5)CO)C6OCC(O)(C6O)CO',
                                     'name': 'Kaempferol '
                                             '3-apiosyl-(1->2)-alpha-L-arabinofuranoside-7-rhamnoside',
                                     'reason': 'No ester group found'},
                                 {   'smiles': '[O-]C(=O)OON=O',
                                     'name': 'nitrosoperoxycarbonate(1-)',
                                     'reason': 'No decanoyl group found'},
                                 {   'smiles': 'O=C1C=2C(=O)OC34C2OC(C(O)C3C=CC(C(CC1C)C)=O)C(C)C4',
                                     'name': 'Atrop-Abybetaomicin C',
                                     'reason': 'No decanoyl group found'},
                                 {   'smiles': 'OC(CC(C)C)/C=C/C',
                                     'name': 'Rhynchophorol.',
                                     'reason': 'No ester group found'},
                                 {   'smiles': 'O=C1C2=C(OC3=C1C=NC(=C3)C)C(OC(=O)C)=C(OC)C=C2C(=O)OC',
                                     'name': 'Diaporphasine C',
                                     'reason': 'No decanoyl group found'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me create a program to classify decanoate esters. The '
                 'key characteristics to check for are:\n'
                 '\n'
                 '1. Presence of a decanoyl group (10-carbon chain with a '
                 'carbonyl)\n'
                 '2. An ester linkage (-C(=O)O-)\n'
                 '3. The decanoyl group must be connected via the ester '
                 'linkage\n'
                 '\n'
                 "Here's the implementation:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 32,
    'num_false_positives': 100,
    'num_true_negatives': 1545,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.24242424242424243,
    'recall': 1.0,
    'f1': 0.3902439024390244,
    'accuracy': 0.9403697078115683,
    'negative_predictive_value': 1.0}