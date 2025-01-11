"""
Classifies: CHEBI:37141 organobromine compound
"""
"""
Classifies: organobromine compound
Definition: A compound containing at least one carbon-bromine bond.
"""
from rdkit import Chem

def is_organobromine_compound(smiles: str):
    """
    Determines if a molecule is an organobromine compound based on its SMILES string.
    An organobromine compound must contain at least one carbon-bromine bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organobromine compound, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of both C and Br atoms
    has_carbon = False
    has_bromine = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon
            has_carbon = True
        elif atom.GetAtomicNum() == 35:  # Bromine
            has_bromine = True
        if has_carbon and has_bromine:
            break
    
    if not has_carbon:
        return False, "No carbon atoms found"
    if not has_bromine:
        return False, "No bromine atoms found"

    # Look for C-Br bond pattern
    cbr_pattern = Chem.MolFromSmarts("[#6]-[Br]")
    if not mol.HasSubstructMatch(cbr_pattern):
        return False, "No carbon-bromine bond found"

    # Count C-Br bonds
    cbr_matches = len(mol.GetSubstructMatches(cbr_pattern))
    
    return True, f"Found {cbr_matches} carbon-bromine bond{'s' if cbr_matches > 1 else ''}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37141',
                          'name': 'organobromine compound',
                          'definition': 'A compound containing at least one '
                                        'carbon-bromine bond.',
                          'parents': ['CHEBI:17792', 'CHEBI:22928'],
                          'xrefs': [   'MetaCyc:Bromide',
                                       'Wikipedia:Organobromine_compound'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [   {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C#CC2=CC=C(C=C2)[C@H]3[C@H](N[C@@H]3C#N)CO)F',
                                     'name': '(2S,3R,4S)-3-[4-[2-(2-fluorophenyl)ethynyl]phenyl]-4-(hydroxymethyl)-2-azetidinecarbonitrile',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'OC(=O)Cc1cn(nc1-c1ccc(Cl)cc1)-c1ccccc1',
                                     'name': 'lonazolac',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O',
                                     'name': 'octadecanedioic acid',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](NC(=O)C)[C@@H]1OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@H]2O)CO[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3S,4R,5R,6R)-5-Acetamido-3,4-dihydroxy-6-[[(2R,3R,4S,5R)-3,4,5,6-tetrahydroxyoxan-2-yl]methoxy]oxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'C\\C=C(/C)CC\\C=C(/C)C(=O)OCC(C)(C)CC1=C(O)C(=O)c2ccccc2C1=O',
                                     'name': 'rhinacanthin C',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'S(=O)(CC1=CC=CC=C1)C',
                                     'name': 'Methyl benzyl sulfoxide',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'C=CCOC1=NS(=O)(=O)c2ccccc12',
                                     'name': 'probenazole',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'C(=O)([C@@H](NC(=O)CC)CCSC)[O-]',
                                     'name': 'N-propanoyl-L-methioninate',
                                     'reason': 'No bromine atoms found'},
                                 {   'smiles': 'O1C(C(O)C(O)C(O)C1OC(CCC=2C(CC(O)CC2C)(C)C)C)COC3OCC(O)(C3O)CO',
                                     'name': '(3b,9R)-5-Megastigmene-3,9-diol '
                                             '9-[apiosyl-(1->6)-glucoside]',
                                     'reason': 'No bromine atoms found'}],
    'sample_false_negatives': [   {   'smiles': 'BrN1C(=O)CCC1=O',
                                      'name': 'N-bromosuccinimide',
                                      'reason': 'No carbon-bromine bond found'},
                                  {   'smiles': 'Brn1ccc2ccccc12',
                                      'name': '1-bromoindole',
                                      'reason': 'No carbon-bromine bond '
                                                'found'}],
    'attempt': 0,
    'reasoning': 'Let me help create a program to classify organobromine '
                 'compounds. The key requirement is to detect at least one '
                 'carbon-bromine bond.\n'
                 '\n'
                 "Here's my step-by-step thought process:\n"
                 '1. First, we need to check if the SMILES string is valid\n'
                 '2. Then, we need to verify that the molecule contains both '
                 'carbon and bromine atoms\n'
                 "3. Finally, we need to check if there's at least one C-Br "
                 'bond\n'
                 '\n'
                 "Here's the code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 308,
    'num_false_positives': 100,
    'num_true_negatives': 12700,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.7549019607843137,
    'recall': 0.9935483870967742,
    'f1': 0.8579387186629526,
    'accuracy': 0.9922196796338673,
    'negative_predictive_value': 0.999842544481184}