"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: gamma-lactone (5-membered lactone ring)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule contains a gamma-lactone (5-membered lactone ring).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a gamma-lactone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Patterns for gamma-lactone detection
    patterns = [
        # Basic gamma-lactone pattern - most general
        # Matches any 5-membered ring with O-C(=O) motif
        "[#8X2r5]1-[#6r5]-[#6r5]-[#6r5]-[#6X3r5](=[O])1",
        
        # Alternative pattern with any bond types
        "[#8X2r5]1~[#6r5]~[#6r5]~[#6r5]~[#6X3r5](=[O])1",
        
        # Pattern for fused systems
        "[#8X2r5](@[#6r5])(@[#6X3r5](=[O]))@[#6r5]@[#6r5]",
        
        # Pattern for substituted variants
        "[#8X2r5]1[#6r5][#6r5][#6r5][#6X3r5](=O)1"
    ]
    
    for pattern in patterns:
        patt = Chem.MolFromSmarts(pattern)
        if patt is not None:
            if mol.HasSubstructMatch(patt):
                # Count unique matches using the most general pattern
                matches = len(mol.GetSubstructMatches(patt))
                if matches == 1:
                    return True, "Contains a gamma-lactone (5-membered lactone ring)"
                else:
                    return True, f"Contains {matches} gamma-lactone rings"
    
    # Additional check for special cases using a more relaxed pattern
    backup_pattern = "[#8r5]1~[#6]~[#6]~[#6]~[#6](=[O])1"
    patt = Chem.MolFromSmarts(backup_pattern)
    if patt is not None and mol.HasSubstructMatch(patt):
        return True, "Contains a gamma-lactone variant"
            
    return False, "No gamma-lactone substructure found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:37581',
                          'name': 'gamma-lactone',
                          'definition': 'A lactone having a five-membered '
                                        'lactone ring.',
                          'parents': ['CHEBI:25000'],
                          'xrefs': ['PMID:18789684'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'CC1CCC2(C(C1(C)CCC3C(C(OC3=O)OC)O)CC(C(C24CO4)O)O)C NAME: '
               "3-[2-(6,7-dihydroxy-1,2,4a-trimethyl-1-spiro[3,4,6,7,8,8a-hexahydro-2H-naphthalene-5,2'-oxirane]yl)ethyl]-4-hydroxy-5-methoxy-2-oxolanone "
               'REASON: MISSED No gamma-lactone substructure found\n'
               ' * SMILES: '
               'O=C1O[C@@H](C)C2=C1[C@H](O)[C@@H](O)[C@@H](O2)/C=C/C NAME: '
               'Mabetaarilactone B REASON: MISSED No gamma-lactone '
               'substructure found\n'
               ' * SMILES: '
               '[H][C@]12CCC3=C[C@@]([H])(OC3=O)c3c(C)coc3C[C@]1(C)O2 NAME: '
               'pseudoneolinderane REASON: MISSED No gamma-lactone '
               'substructure found\n'
               ' * SMILES: '
               'OC12C3C(C4(C(O)(CC3)CC(OC5OC(C(OC6OC(C(OC7OC(C(O)C(O)C7O)CO)C(O)C6O)CO)C(O)C5)C)CC4)CO)CCC1(C(CC2)C=8COC(=O)C8)C '
               'NAME: Corchorusoside E REASON: MISSED No gamma-lactone '
               'substructure found\n'
               ' * SMILES: O1C(CCC1=O)CC2=CC=CC=C2 NAME: 5-benzyloxolan-2-one '
               'REASON: MISSED No gamma-lactone substructure found\n'
               ' * SMILES: '
               'O=C1OC(=O)C(=C1/C=C/CC)CC(CC=2C(=O)OC(C2CC(CC=3C(=O)OC(C3CCC(=O)O)=O)CC)=O)CC '
               'NAME: Cordyanhydride B REASON: MISSED No gamma-lactone '
               'substructure found\n'
               ' * SMILES: O=C1O[C@H](CCCCCCCCCCCCCC(=O)C)C(C1=C)C(=O)O NAME: '
               'Muronic acid REASON: MISSED No gamma-lactone substructure '
               'found\n'
               ' * SMILES: O=C1OC(CCCCCCCCCCCCCC(=O)C)C(C1=C)C(=O)O NAME: '
               'Pertusaric acid REASON: MISSED No gamma-lactone substructure '
               'found\n'
               ' * SMILES: CC(SC#N)C1CC(C)(C)OC1=O NAME: '
               '5,5-dimethyl-3-(1-thiocyanatoethyl)dihydro-2(3H)-furanone '
               'REASON: MISSED No gamma-lactone substructure found\n'
               ' * SMILES: CCCCC1COC(=O)C1 NAME: 4-butyl-gamma-butyrolactone '
               'REASON: MISSED No gamma-lactone substructure found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C#CC2=CC=C(C=C2)[C@H]3[C@H](N[C@@H]3C#N)CO)F',
                                     'name': '(2S,3R,4S)-3-[4-[2-(2-fluorophenyl)ethynyl]phenyl]-4-(hydroxymethyl)-2-azetidinecarbonitrile',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'OC(=O)Cc1cn(nc1-c1ccc(Cl)cc1)-c1ccccc1',
                                     'name': 'lonazolac',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O',
                                     'name': 'octadecanedioic acid',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'O1[C@@H]([C@@H](O)[C@H](O)[C@@H](NC(=O)C)[C@@H]1OC[C@H]2OC(O)[C@H](O)[C@@H](O)[C@H]2O)CO[C@@H]3O[C@@H]([C@@H](O)[C@H](O)[C@H]3NC(=O)C)CO',
                                     'name': 'N-[(2R,3R,4R,5S,6R)-2-[[(2R,3S,4R,5R,6R)-5-Acetamido-3,4-dihydroxy-6-[[(2R,3R,4S,5R)-3,4,5,6-tetrahydroxyoxan-2-yl]methoxy]oxan-2-yl]methoxy]-4,5-dihydroxy-6-(hydroxymethyl)oxan-3-yl]acetamide',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'C\\C=C(/C)CC\\C=C(/C)C(=O)OCC(C)(C)CC1=C(O)C(=O)c2ccccc2C1=O',
                                     'name': 'rhinacanthin C',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C(=C3)Br)O)Br)OC',
                                     'name': '2,6-dibromo-4-[(6,7-dimethoxy-4-quinazolinyl)amino]phenol',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'S(=O)(CC1=CC=CC=C1)C',
                                     'name': 'Methyl benzyl sulfoxide',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'C=CCOC1=NS(=O)(=O)c2ccccc12',
                                     'name': 'probenazole',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'},
                                 {   'smiles': 'C(=O)([C@@H](NC(=O)CC)CCSC)[O-]',
                                     'name': 'N-propanoyl-L-methioninate',
                                     'reason': 'No gamma-lactone substructure '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': '[H]C(OCC)=C1N=C(OC1=O)c1ccco1',
                                      'name': '4-(ethoxymethylene)-2-(2-furyl)oxazol-5-one',
                                      'reason': 'No gamma-lactone substructure '
                                                'found'},
                                  {   'smiles': '[H]C(O)=C1N=C(Cc2ccccc2)OC1=O',
                                      'name': '2-benzyl-4-hydroxymethylene-5-oxazolone',
                                      'reason': 'No gamma-lactone substructure '
                                                'found'},
                                  {   'smiles': '[H]C(OCC)=C1N=C(OC1=O)c1ccccc1',
                                      'name': '4-(ethoxymethylene)-2-phenyloxazol-5-one',
                                      'reason': 'No gamma-lactone substructure '
                                                'found'},
                                  {   'smiles': '[H]C(OCC)=C1N=C(OC1=O)c1ccc(cc1)[N+]([O-])=O',
                                      'name': '2-(p-nitrophenyl)-4-ethoxymethyleneoxazol-5-one',
                                      'reason': 'No gamma-lactone substructure '
                                                'found'},
                                  {   'smiles': '[H][C@]12CC[C@](C)(O)[C@@H](CC3=CC(=O)C=C(OC)C3=O)[C@]1(C)CC[C@@]1([H])[C@]3(C)CCC[C@]21COC3=O',
                                      'name': 'strongylophorine-26',
                                      'reason': 'No gamma-lactone substructure '
                                                'found'},
                                  {   'smiles': '[H]C(OCC)=C1N=C(OC1=O)\\C=C\\c1ccccc1',
                                      'name': '4-(ethoxymethylene)-2-styryloxazol-5-one',
                                      'reason': 'No gamma-lactone substructure '
                                                'found'},
                                  {   'smiles': 'C[C@@H](O)[C@@H]1N\\C(OC1=O)=N/c1ncnc2n(cnc12)[C@@H]1O[C@H](CO)[C@@H](O)[C@H]1O',
                                      'name': 'cyclic '
                                              'N(6)-threonylcarbamoyladenosine',
                                      'reason': 'No gamma-lactone substructure '
                                                'found'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous code and the failed outcomes:\n'
                 '\n'
                 '1. The main issue appears to be that the SMARTS patterns '
                 'were too restrictive. The patterns were looking for very '
                 'specific arrangements of atoms, but gamma-lactones can have '
                 'various substituents and configurations.\n'
                 '\n'
                 '2. Looking at the failed examples:\n'
                 '- Many contain gamma-lactone rings with additional '
                 'substituents or fused ring systems\n'
                 '- Some have stereochemistry which might have interfered with '
                 'matching\n'
                 "- The patterns didn't account for variations in bond types "
                 'within the ring\n'
                 '\n'
                 '3. Key improvements needed:\n'
                 '- Make patterns more flexible to match substituted '
                 'gamma-lactones\n'
                 '- Better handle fused ring systems\n'
                 '- Ignore stereochemistry in the initial matching\n'
                 '- Consider using simpler patterns with fewer constraints\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 966,
    'num_false_positives': 100,
    'num_true_negatives': 4535,
    'num_false_negatives': 7,
    'num_negatives': None,
    'precision': 0.9061913696060038,
    'recall': 0.9928057553956835,
    'f1': 0.9475232957332026,
    'accuracy': 0.9809201141226819,
    'negative_predictive_value': 0.9984588287098195}