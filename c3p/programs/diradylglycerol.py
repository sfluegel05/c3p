"""
Classifies: CHEBI:76578 diradylglycerol
"""
from rdkit import Chem

def is_diradylglycerol(smiles: str):
    """
    Determines if a molecule is a diradylglycerol based on its SMILES string.
    A diradylglycerol is characterized by a glycerol backbone with two fatty acid chains attached as esters or ethers,
    and one free hydroxyl group (usually at the sn-3 position).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Correct glycerol backbone pattern C(CO)CO with potential stereochemistry
    glycerol_pattern = Chem.MolFromSmarts("OCC(O)CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Identify ester: R-C(=O)-O or ether: R-O-R' linkages
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ether_pattern = Chem.MolFromSmarts("COC")
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)

    # Total should be at least 2 esters/ethers indicating two fatty acid chains
    if len(ester_matches) + len(ether_matches) < 2:
        return False, f"Found {len(ester_matches)} esters and {len(ether_matches)} ethers, need at least two in any combination"
    
    return True, "Contains glycerol backbone with two ester or ether linkages"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:76578',
                          'name': 'diradylglycerol',
                          'definition': 'Any lipid that is glycerol bearing '
                                        'two substituent groups - either acyl, '
                                        'alkyl, or alk-1-enyl - at any two of '
                                        'the three possible positions.',
                          'parents': ['CHEBI:35741'],
                          'xrefs': ['LIPID_MAPS_class:LMGL0201'],
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
    'sample_true_negatives': [   {   'smiles': 'CC1=CC(=C(C=C1)C)C(=O)CSC2=NN=C(S2)C',
                                     'name': '1-(2,5-dimethylphenyl)-2-[(5-methyl-1,3,4-thiadiazol-2-yl)thio]ethanone',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C1C=C([C@@]2(C(C[C@@H](C2)O)(C)C)C)CC[C@@]1(O)C',
                                     'name': 'Enokipodin H',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'CCCC[C@](Cn1cncn1)(C#N)c1ccc(Cl)cc1',
                                     'name': '(S)-myclobutanil',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O=C1O[C@@H]([C@H](NC(=O)[C@@H](NC(=O)CCCCC)CC(=O)O)C(=O)N[C@H](C(=O)N[C@H]2CC[C@H](N([C@H](C(N([C@H](C(N[C@H]1C(C)C)=O)CC3=CC=CC=C3)C)=O)CC(C)C)C2=O)O)CCCCN(C)C)C',
                                     'name': 'Cyanopeptolin D',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'O1[C@]2(O)N([C@H](C(=O)N3[C@]2(CCC3)[H])C(C)C)C(=O)[C@@]1(NC(=O)[C@H]4CN([C@]5(C(=C4)C6=C7C(C5)=CNC7=CC=C6)[H])C)C',
                                     'name': 'Ergovaline',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'OC(=O)CCCCCCC#CCCCC',
                                     'name': '8-tridecynoic acid',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'CC(=O)N([O-])CCC[C@H](N)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CCCN([O-])C(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@H]([C@H](O)[C@H]1S[C@H]([C@H](O)[C@H]1O)n1ccc(=N)n(C)c1=O)C(O)=O',
                                     'name': 'desferrialbomycin epsilon(3-)',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'SC[C@H](N)C(=O)N[C@H](C(=O)N[C@@H](CCC(=O)N)C(O)=O)CO',
                                     'name': 'Cys-Ser-Gln',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'S(C(C)C(O)=O)CC1=CC=CC=C1',
                                     'name': '2-(Benzylthio)propanoic acid',
                                     'reason': 'No glycerol backbone found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)NC3=CC=C(C=C3)Cl)C(=O)N(C[C@@H]1OC)C)C)C',
                                     'name': '1-(4-chlorophenyl)-3-[(4S,7R,8R)-8-methoxy-4,5,7,10-tetramethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'No glycerol backbone found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': 'Upon reviewing the previous program and the outcomes, the '
                 'primary issue seems to be with the definition and matching '
                 'of the glycerol backbone pattern. The SMARTS pattern '
                 '"[C@@H](CO)CO | [C@H](CO)CO" was used to identify the '
                 'glycerol backbone, but it appears that this pattern is '
                 'incorrect and a potential source of the error, as no true '
                 'positives were identified.\n'
                 '\n'
                 'To solve this issue, we need to accurately define the '
                 'criteria that characterize a diradylglycerol molecule. Key '
                 'defining features for diradylglycerols include:\n'
                 '- A glycerol backbone (1,2,3-propanetriol) where only two '
                 'hydroxyl groups are esterified or etherified with fatty '
                 'acids, leaving one free hydroxyl group.\n'
                 '- The ester or ether linkages to fatty acid chains typically '
                 'found at positions 1 and 2 or 1 and 3.\n'
                 '\n'
                 '1. We start by verifying the glycerol core (C-C-C with '
                 'hydroxyl groups), usually represented as C[C@@H](O)CO or '
                 'C[C@H](O)CO, while accounting for the possible '
                 'stereochemistry at the central carbon.\n'
                 '\n'
                 '2. Next, we look for the presence of ester and/or ether '
                 'linkages. Ester linkage can be identified by the SMILES '
                 'pattern "C(=O)O", and ethers by the pattern "COC".\n'
                 '\n'
                 '3. We check if there are two ester or ether linkages '
                 'combined, identifying them on either the 1,2 or 1,3 '
                 'positions of the glycerol backbone.\n'
                 '\n'
                 "Let's correct the code by updating the SMARTS patterns and "
                 'logic:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 677,
    'num_false_positives': 100,
    'num_true_negatives': 442,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.8712998712998713,
    'recall': 1.0,
    'f1': 0.9312242090784043,
    'accuracy': 0.9179655455291222,
    'negative_predictive_value': 1.0}