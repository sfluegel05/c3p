"""
Classifies: CHEBI:24128 furanocoumarin
"""
from rdkit import Chem

def is_furanocoumarin(smiles: str):
    """
    Determines if a molecule is a furanocoumarin based on its SMILES string.
    A furanocoumarin is defined as a furochromene that consists of a furan ring 
    fused with a coumarin.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a furanocoumarin, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string to a molecular object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for a furanocoumarin backbone
    # Patterns should capture the various potential fusion types of furo and coumarin
    furanocoumarin_patterns = [
        Chem.MolFromSmarts("c1oc2ccccc2oc1C=O"), # Basic furo[2,3-b]coumarin structure
        Chem.MolFromSmarts("c1oc2cccc(c2cc1C=O)"), # Basic furo[3,2-g]coumarin structure
        Chem.MolFromSmarts("c1oc2ccc(o2)c(c1)C=O")  # Basic psoralen-like structure
    ]
    
    # Ensure all patterns are valid
    if any(pattern is None for pattern in furanocoumarin_patterns):
        return None, None  # If any pattern is None, there's a problem with its definition
    
    # Check for presence of any furanocoumarin patterns
    for pattern in furanocoumarin_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains furanocoumarin fused rings"
            
    return False, "The structure lacks necessary furanocoumarin fused moieties"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:24128',
                          'name': 'furanocoumarin',
                          'definition': 'Any furochromene that consists of a '
                                        'furan ring fused with a coumarin. The '
                                        'fusion may occur in different ways in '
                                        'give several isomers.',
                          'parents': ['CHEBI:39432'],
                          'xrefs': ['Wikipedia:Furocoumarin'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Error: Python argument types in\n'
               '    Mol.HasSubstructMatch(Mol, NoneType)\n'
               'did not match C++ signature:\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, RDKit::SubstructMatchParameters params=True)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'RDKit::SubstructMatchParameters params)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::MolBundle '
               'query, bool recursionPossible=True, bool useChirality=False, '
               'bool useQueryQueryMatches=False)\n'
               '    HasSubstructMatch(RDKit::ROMol self, RDKit::ROMol query, '
               'bool recursionPossible=True, bool useChirality=False, bool '
               'useQueryQueryMatches=False)\n'
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
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'O1C2=C(C=3OC(=O)C=CC3C=C2)C=C1C(=O)C',
                                      'name': "2'-Acetylangelicin",
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'O1C=2C(=C(OC)C3=C(OC(=O)C=C3)C2O)C=C1',
                                      'name': '9-Hydroxy-4-methoxypsoralen',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'CC(C)=CC(O)C\\C(C)=C\\COc1c2ccoc2cc2oc(=O)ccc12',
                                      'name': 'Notopterol',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'O1C(C(OC(=O)/C=C/C2=CC=C(O)C=C2)(C)C)C(O)C3=C1C=CC4=C3OC(=O)C=C4',
                                      'name': '2-{9-hydroxy-2-oxo-2H,8H,9H-furo[2,3-h]chromen-8-yl}propan-2-yl '
                                              '(2E)-3-(4-hydroxyphenyl)prop-2-enoate',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'O1C(CC=2C1=C3C(OC(=O)C=C3CCC)=C(C2O)C(=O)C(CC)C)C(O)(C)C',
                                      'name': 'Cycloneomammein',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'COc1c2ccoc2cc2oc(=O)ccc12',
                                      'name': '5-methoxypsoralen',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'CC[C@H](C)C(=O)c1c(O)c2c(cc(=O)oc2c2cc(oc12)C(C)(C)O)-c1ccccc1',
                                      'name': 'ochrocarpin A',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'CC(C)(O)[C@H](O)COc1c2ccoc2cc2oc(=O)ccc12',
                                      'name': 'oxypeucedanin hydrate',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'Cc1cc2occc2c2oc(=O)ccc12',
                                      'name': '5-Methylangelicin',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'},
                                  {   'smiles': 'CC(=O)O[C@H]1[C@H](Oc2ccc3ccc(=O)oc3c12)C(C)(C)OC(=O)C=C(C)C',
                                      'name': 'Peucenidin',
                                      'reason': 'The structure lacks necessary '
                                                'furanocoumarin fused '
                                                'moieties'}],
    'attempt': 2,
    'reasoning': 'In order to create a program that accurately classifies '
                 'furanocoumarins based on their SMILES strings, we should '
                 "first address the previous code's issues. The error message "
                 'suggests that the substructure match failed due to the '
                 'pattern being `None`. This likely means the SMARTS patterns '
                 'were invalid or not appropriate for capturing the structures '
                 'we are targeting.\n'
                 '\n'
                 'Furanocoumarins are characterized by a fusion of a furan '
                 'ring with a coumarin, and they have variability in their '
                 "fusion (different isomers). Let's break down our approach:\n"
                 '\n'
                 '1. **Review SMARTS Patterns**: We need detailed SMARTS '
                 'patterns to capture the fusion of furan and coumarin '
                 'correctly. The patterns in the previous attempt may have '
                 'been too simplistic or incorrect.\n'
                 '\n'
                 '2. **Fusion Sites and Isomers**: We should ensure that the '
                 'SMARTS captures different possible ways a furan can be fused '
                 'with a coumarin, accounting for isomeric variation.\n'
                 '\n'
                 '3. **Validation and Testing**: Use the example SMILES '
                 'strings provided to test that our function can identify a '
                 'variety of furanocoumarins.\n'
                 '\n'
                 "Let's refine the code to improve recognition of "
                 'furanocoumarins:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 2,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.0196078431372549,
    'f1': 0.038461538461538464,
    'accuracy': 0.0196078431372549,
    'negative_predictive_value': 0.0}