"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:52247 octanoate ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    An octanoate ester is any ester where the acid component is octanoic acid (caprylic acid).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group pattern
    ester_pattern = Chem.MolFromSmarts("[#6][CX3](=[OX1])[OX2][#6]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"
    
    # Pattern for octanoate group: linear 7-carbon chain attached to ester carbonyl
    # Note: The SMARTS below looks for:
    # - An ester group (-C(=O)O-)
    # - Connected to exactly 6 more carbons in a linear chain
    # - The last carbon must have 3 hydrogens (CH3)
    octanoate_pattern = Chem.MolFromSmarts("[CH3][CH2][CH2][CH2][CH2][CH2][CH2]C(=O)[OX2]")
    
    if not mol.HasSubstructMatch(octanoate_pattern):
        return False, "No octanoate (8-carbon) chain found attached to ester group"
    
    # Count matches to ensure we have at least one octanoate ester
    matches = mol.GetSubstructMatches(octanoate_pattern)
    
    # Get the atoms involved in each match to verify the chain length
    for match in matches:
        # The pattern above ensures we have the right structure
        return True, "Contains octanoate ester group (8-carbon chain with terminal ester)"
        
    return False, "No valid octanoate ester found"

def test_examples():
    """Test function with some known examples"""
    test_cases = [
        ("CCCCCCCC(=O)OCC", True),  # ethyl octanoate
        ("CCCCCCCC(=O)OC", True),   # methyl octanoate
        ("CCCCCCC(=O)OCC", False),  # ethyl heptanoate (too short)
        ("CCCCCCCCC(=O)OCC", False),# ethyl nonanoate (too long)
        ("CCCCCCCCOC", False),      # octyl ether (not an ester)
    ]
    
    for smiles, expected in test_cases:
        result, reason = is_octanoate_ester(smiles)
        print(f"SMILES: {smiles}")
        print(f"Expected: {expected}, Got: {result}")
        print(f"Reason: {reason}\n")

if __name__ == "__main__":
    test_examples()


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:87657',
                          'name': 'octanoate ester',
                          'definition': 'Any fatty acid ester in which the '
                                        'carboxylic acid component is octanoic '
                                        'acid (caprylic acid).',
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
                                     'reason': 'No ester group found'},
                                 {   'smiles': 'P(OCC(O)COC(=O)CCCCCCCCCCCCCCCCC)(O)(O)=O',
                                     'name': 'LysoPA(18:0/0:0)',
                                     'reason': 'No octanoate (8-carbon) chain '
                                               'found attached to ester group'},
                                 {   'smiles': 'O=C1C=2C(=O)OC34C2OC(C(O)C3C=CC(C(CC1C)C)=O)C(C)C4',
                                     'name': 'Atrop-Abybetaomicin C',
                                     'reason': 'No octanoate (8-carbon) chain '
                                               'found attached to ester group'},
                                 {   'smiles': 'P(OC[C@H](OC(=O)CCCCCCCCCCCCCCCCCC)CO/C=C\\CCCCCCCCCCCCCCCC)(OC[C@@H](O)CO)(O)=O',
                                     'name': 'PG(P-18:0/19:0)',
                                     'reason': 'No octanoate (8-carbon) chain '
                                               'found attached to ester '
                                               'group'}],
    'sample_false_negatives': [],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify octanoate esters. '
                 'The key characteristics are:\n'
                 '\n'
                 '1. Must contain an ester group (-C(=O)O-)\n'
                 '2. One part of the ester must be octanoic acid (C8H16O2, or '
                 'CH3(CH2)6COOH before esterification)\n'
                 '3. The other part can be any group (alcohol part)\n'
                 '\n'
                 "Here's the program:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 28,
    'num_false_positives': 100,
    'num_true_negatives': 81489,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.21875,
    'recall': 1.0,
    'f1': 0.358974358974359,
    'accuracy': 0.9987747650612004,
    'negative_predictive_value': 1.0}