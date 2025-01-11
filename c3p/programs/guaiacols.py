"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: guaiacols
Definition: Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_guaiacols(smiles: str):
    """
    Determines if a molecule contains a guaiacol moiety based on its SMILES string.
    A guaiacol is a phenol with a methoxy group in the ortho position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains guaiacol moiety, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # More flexible SMARTS patterns for guaiacol core:
    # Match aromatic carbon with OH and adjacent carbon with OMe
    # Try both orientations (OH-OMe and OMe-OH)
    patterns = [
        # Pattern 1: OH then OMe
        '[c;$(c1ccccc1)]([OH1])[$(c([OC]))]-,=',
        # Pattern 2: OMe then OH
        '[c;$(c1ccccc1)]([OCH3])[$(c([OH1]))]-,='
    ]
    
    for pattern in patterns:
        guaiacol_pattern = Chem.MolFromSmarts(pattern)
        if guaiacol_pattern is None:
            continue
            
        matches = mol.GetSubstructMatches(guaiacol_pattern)
        if matches:
            # Verify each match
            for match in matches:
                # Get the matched atoms
                core_atoms = [mol.GetAtomWithIdx(idx) for idx in match]
                
                # Check if the carbons are part of an aromatic ring
                ring_info = mol.GetRingInfo()
                if all(ring_info.IsAtomInRingOfSize(atom.GetIdx(), 6) for atom in core_atoms 
                      if atom.GetAtomicNum() == 6):
                    return True, "Contains phenol with ortho methoxy group (guaiacol moiety)"
    
    # Additional check with a simpler pattern
    simple_pattern = Chem.MolFromSmarts('c(O)c(OC)')
    if mol.HasSubstructMatch(simple_pattern):
        return True, "Contains phenol with ortho methoxy group (guaiacol moiety)"
                
    return False, "No guaiacol pattern found"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:134251',
                          'name': 'guaiacols',
                          'definition': 'Any phenol carrying an additional '
                                        'methoxy substituent at the '
                                        'ortho-position.',
                          'parents': ['CHEBI:25235', 'CHEBI:33853'],
                          'xrefs': ['MetaCyc:Guaiacols'],
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
               'C=1C(=C(C=CC1/C=C/C(O)=O)OC(CO)C(OS(O)(=O)=O)C=2C=C(C(=CC2)O)OC)OC '
               'NAME: '
               '(2E)-3-(4-{[3-hydroxy-1-(4-hydroxy-3-methoxyphenyl)-1-(sulfooxy)propan-2-yl]oxy}-3-methoxyphenyl)prop-2-enoic '
               'acid REASON: MISSED Missing required elements (C, O)\n'
               ' * SMILES: COc1cc(\\C=C/C(=O)NCCCCNC(N)=N)ccc1O NAME: '
               '(Z)-4-hydroxy-3-methoxycinnamoylagmatine REASON: MISSED '
               'Missing required elements (C, O)\n'
               ' * SMILES: COc1cc(\\C=C\\COC(=O)\\C=C\\c2ccc(O)cc2)ccc1O NAME: '
               'coniferyl p-coumarate REASON: MISSED Contains required '
               'elements but not in correct configuration\n'
               ' * SMILES: COc1cc(CC(O)=O)ccc1O NAME: homovanillic acid '
               'REASON: MISSED Contains required elements but not in correct '
               'configuration\n'
               ' * SMILES: '
               'C1=C(C=C2C(=C1OC)O[C@H]([C@@H]2C(=O)O)C3=CC=C(C(=C3)OC)O)/C=C/C(O)=O '
               'NAME: (-)-DCA-CC REASON: MISSED Contains required elements but '
               'not in correct configuration\n'
               ' * SMILES: '
               'OC1=C(C\\C=C(\\CC\\C=C(\\CC\\C=C(\\CC\\C=C(\\CC\\C=C(\\CC\\C=C(\\CC\\C=C(\\CCC=C(C)C)/C)/C)/C)/C)/C)/C)/C)C=CC=C1OC '
               'NAME: 2-methoxy-6-(all-trans-octaprenyl)phenol REASON: MISSED '
               'Contains required elements but not in correct configuration\n'
               ' * SMILES: COc1cc(cc(O)c1O)C(=O)Cc1cccc(O)c1C(O)=O NAME: '
               '2-hydroxy-6-[2-(3,4-dihydroxyphenyl-5-methoxy)-2-oxoethyl]benzoic '
               'acid REASON: MISSED Contains required elements but not in '
               'correct configuration\n'
               ' * SMILES: '
               'COc1cc(ccc1O)[C@@H]1[C@@H](CO)[C@H](CO)Cc2cc(OC)c(O)cc12 NAME: '
               '(+)-isolariciresinol REASON: MISSED Contains required elements '
               'but not in correct configuration\n'
               ' * SMILES: '
               'COc1cc(C[C@H]2CO[C@@H]([C@H]2CO)c2ccc(O)c(O)c2)ccc1O NAME: '
               'taxiresinol REASON: MISSED Contains required elements but not '
               'in correct configuration\n'
               ' * SMILES: CCCCc1ccc(CCNC(=O)\\C=C\\c2ccc(O)c(OC)c2)cc1 NAME: '
               'Cuscuta propenamide 2 REASON: MISSED Missing required elements '
               '(C, O)\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'N12CN3CN(CN(C1)CC3)CC2',
                                     'name': '1,3,6,8-tetraazatricyclo[4,4,1,1(3,8)]dodecane',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': 'CC1CC2C(CC=C1CCO)C(C)C(=O)O2',
                                     'name': '(8alpha,10beta,11beta)-3-hydroxy-4,15-dinor-1(5)-xanthen-12,8-olide',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(O)=O)[C@H](O)\\C=C\\CCCCCCCCCCCCC',
                                     'name': 'N-docosanoylsphingosine-1-phosphate',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': 'C1CNCCC1SC2=NC(=CC=C2)Cl',
                                     'name': '2-chloro-6-(4-piperidinylthio)pyridine',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': '[H][C@@]1(CC[C@@]2(C)[C@]3([H])CC[C@@]4([H])[C@@](C)([C@H](CC[C@@]44C[C@@]34CC[C@]12C)OC(C)=O)C(O)=O)[C@H](C)CCC(=O)C(C)C',
                                     'name': 'bonianic acid B, (rel)-',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': 'O=C1OC[C@@]2(C1=C(CC[C@]34[C@H]2C[C@H](CC[C@H]3C)C4(C)C)C)C',
                                     'name': 'Harzianelactone',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': 'O=C(NC1=NC=C(Cl)C=C1)CSC2=NN=C3C(NC4=C3C=CC=C4)=N2',
                                     'name': 'dCeMM2',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': 'N[C@@H](CS[C@H](\\C=C\\C=C\\C=C/C=C/CC(O)=O)[C@@H](O)CCCC(O)=O)C(O)=O',
                                     'name': '(13E)-16-carboxy-Delta(13)-17,18,19,20-tetranor-leukotriene '
                                             'E4',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](OC[C@H]2OC(O)[C@H](NC(C)=O)[C@@H](O[C@@H]3O[C@@H](C)[C@@H](O)[C@@H](O)[C@@H]3O)[C@@H]2O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->3)-[alpha-L-Fucp-(1->6)]-{beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-[beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->6)]-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)}-D-GlcpNAc',
                                     'reason': 'No guaiacol pattern found'},
                                 {   'smiles': 'CCCCCCCCCCCCCCCCCC(=O)OC[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1ccc(N)nc1=O)OC(=O)CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC',
                                     'name': 'CDP-1-stearoyl-2-(4Z,7Z,10Z,13Z,16Z,19Z)-docosahexaenoyl-sn-glycerol',
                                     'reason': 'No guaiacol pattern found'}],
    'sample_false_negatives': [],
    'attempt': 3,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. The main issue seems to be that the SMARTS pattern was '
                 'too restrictive. The pattern `c1(O)c(OC)cccc1` requires a '
                 'complete benzene ring with specific positions for the OH and '
                 'OMe groups.\n'
                 '\n'
                 '2. Many of the positive examples have substitutions on the '
                 'benzene ring, which might interfere with the matching of the '
                 'complete ring pattern.\n'
                 '\n'
                 '3. The code was also checking for aromaticity too strictly - '
                 'while guaiacols are aromatic, RDKit might not always encode '
                 'this correctly in complex molecules.\n'
                 '\n'
                 "4. The element checking was unnecessary since we're "
                 'specifically looking for phenol and methoxy groups which '
                 'must contain C and O.\n'
                 '\n'
                 "Here's an improved version that should be more robust:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 81,
    'num_false_positives': 100,
    'num_true_negatives': 1850,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.44751381215469616,
    'recall': 1.0,
    'f1': 0.6183206106870229,
    'accuracy': 0.9507631708517972,
    'negative_predictive_value': 1.0}