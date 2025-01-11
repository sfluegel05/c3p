"""
Classifies: CHEBI:29017 1-acyl-sn-glycero-3-phosphoethanolamine
"""
"""
Classifies: 1-acyl-sn-glycero-3-phosphoethanolamine
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_1_acyl_sn_glycero_3_phosphoethanolamine(smiles: str):
    """
    Determines if a molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 1-acyl-sn-glycero-3-phosphoethanolamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for basic elements (excluding H since it might be implicit)
    required_elements = {'C', 'O', 'N', 'P'}
    mol_elements = set(atom.GetSymbol() for atom in mol.GetAtoms())
    if not required_elements.issubset(mol_elements):
        return False, f"Missing required elements. Found {mol_elements}, need {required_elements}"

    # Check for glycerol backbone with correct stereochemistry
    # Multiple SMARTS patterns to catch different representations
    glycerol_patterns = [
        "[CH2X4][C@@HX4][CH2X4]",  # R configuration
        "[CH2X4][C@HX4][CH2X4]",   # Alternative representation
        "[CH2][C@H][CH2]",         # Simplified pattern
        "[CH2][C@@H][CH2]"         # Simplified alternative
    ]
    
    found_glycerol = False
    for pattern in glycerol_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_glycerol = True
            break
    
    if not found_glycerol:
        return False, "No glycerol backbone with correct stereochemistry found"

    # Check for one ester group (acyl chain at position 1)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])[#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check for phosphoethanolamine group - multiple patterns to catch different forms
    phosphoethanolamine_patterns = [
        "[OX2][PX4](=[OX1])([OX2])[OX2]CC[NX3]",           # Neutral form
        "[OX2][PX4](=[OX1])([OX2])[OX2]CC[NH3+]",         # Protonated amine
        "[OX2][PX4](=[OX1])([O-])[OX2]CC[NH3+]",          # Zwitterionic form
        "[OX2][PX4](=[OX1])([OX2H])[OX2]CC[NX3]"          # With explicit H
    ]
    
    found_pe = False
    for pattern in phosphoethanolamine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_pe = True
            break
            
    if not found_pe:
        return False, "No phosphoethanolamine group found"

    # Check for hydroxyl group - multiple patterns
    hydroxyl_patterns = [
        "[OX2H1]",            # Explicit H
        "[OX2H]",             # Alternative representation
        "[OH]"                # Simplified
    ]
    
    found_hydroxyl = False
    for pattern in hydroxyl_patterns:
        matches = mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))
        if len(matches) >= 1:
            found_hydroxyl = True
            break
            
    if not found_hydroxyl:
        return False, "No free hydroxyl group found"

    # Check carbon chain length (should be at least 13 carbons total including glycerol backbone)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 13:
        return False, f"Carbon count too low ({c_count}), need at least 13"

    return True, "Contains glycerol backbone with correct stereochemistry, one acyl chain, phosphoethanolamine group, and free hydroxyl group"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:29017',
                          'name': '1-acyl-sn-glycero-3-phosphoethanolamine',
                          'definition': 'A 1-O-acylglycerophosphoethanolamine '
                                        'having (R)-configuration.',
                          'parents': ['CHEBI:55493'],
                          'xrefs': [   'KEGG:C04438',
                                       'LIPID_MAPS_instance:LMGP02050000'],
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
               'CCCCCCCC\\C=C/CCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN NAME: '
               '1-oleoyl-sn-glycero-3-phosphoethanolamine REASON: MISSED '
               "Missing required elements. Found {'P', 'C', 'O', 'N'}, need "
               "{'H', 'P', 'C', 'O', 'N'}\n"
               ' * SMILES: [C@@H](COC(=O)CCCCCCCCCCCCCCC)(COP(OCCN)(=O)O)O '
               'NAME: 1-hexadecanoyl-sn-glycero-3-phosphoethanolamine REASON: '
               "MISSED Missing required elements. Found {'P', 'C', 'O', 'N'}, "
               "need {'H', 'P', 'C', 'O', 'N'}\n"
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CC/C=C\\C/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(OCCN)(O)=O '
               'NAME: LysoPE(22:5(4Z,7Z,10Z,13Z,16Z)/0:0) REASON: MISSED '
               "Missing required elements. Found {'C', 'P', 'O', 'N'}, need "
               "{'H', 'P', 'C', 'O', 'N'}\n"
               ' * SMILES: '
               '[C@](COC(=O)CCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CCCCC)(O)([H])COP(OCCN)(O)=O '
               'NAME: PE(22:4(7Z,10Z,13Z,16Z)/0:0) REASON: MISSED Missing '
               "required elements. Found {'P', 'C', 'O', 'N'}, need {'H', 'P', "
               "'C', 'O', 'N'}\n"
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN '
               'NAME: 1-arachidonoyl-sn-glycero-3-phosphoethanolamine REASON: '
               "MISSED Missing required elements. Found {'P', 'C', 'O', 'N'}, "
               "need {'H', 'P', 'C', 'O', 'N'}\n"
               ' * SMILES: P(OC[C@H](O)COC(=O)CCCCCCCCCCCC)(OCCN)(O)=O NAME: '
               'PE(13:0/0:0) REASON: MISSED Missing required elements. Found '
               "{'C', 'P', 'O', 'N'}, need {'H', 'P', 'C', 'O', 'N'}\n"
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCCCC/C=C\\CCCCCCCCCC)(OCCN)(O)=O NAME: '
               'PE(22:1(11Z)/0:0) REASON: MISSED Missing required elements. '
               "Found {'C', 'P', 'O', 'N'}, need {'H', 'P', 'C', 'O', 'N'}\n"
               ' * SMILES: '
               'P(OC[C@H](O)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)(OCCN)(O)=O '
               'NAME: LysoPE(18:3(9Z,12Z,15Z)/0:0) REASON: MISSED Missing '
               "required elements. Found {'C', 'P', 'O', 'N'}, need {'H', 'P', "
               "'C', 'O', 'N'}\n"
               ' * SMILES: C(CN)OP(=O)(O)OC[C@H](O)COC(CCCCCCCCCCCCCCCCCCC)=O '
               'NAME: 1-icosanoyl-sn-glycero-3-phosphoethanolamine REASON: '
               "MISSED Missing required elements. Found {'P', 'C', 'O', 'N'}, "
               "need {'H', 'P', 'C', 'O', 'N'}\n"
               ' * SMILES: CCCCCCCCCCCCCCCCC(=O)OC[C@@H](O)COP(O)(=O)OCCN '
               'NAME: 1-heptadecanoyl-sn-glycero-3-phosphoethanolamine REASON: '
               "MISSED Missing required elements. Found {'P', 'C', 'O', 'N'}, "
               "need {'H', 'P', 'C', 'O', 'N'}\n"
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'Missing required elements. '
                                               "Found {'C'}, need {'P', 'C', "
                                               "'O', 'N'}"},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'Missing required elements. '
                                               "Found {'C', 'O'}, need {'P', "
                                               "'C', 'O', 'N'}"},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Missing required elements. '
                                               "Found {'C', 'O', 'N'}, need "
                                               "{'P', 'C', 'O', 'N'}"},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'Missing required elements. '
                                               "Found {'C', 'O', 'N'}, need "
                                               "{'P', 'C', 'O', 'N'}"},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'Missing required elements. '
                                               "Found {'C', 'O', 'N'}, need "
                                               "{'P', 'C', 'O', 'N'}"},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'Missing required elements. '
                                               "Found {'C', 'O'}, need {'P', "
                                               "'C', 'O', 'N'}"},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'Missing required elements. '
                                               "Found {'C', 'O'}, need {'P', "
                                               "'C', 'O', 'N'}"},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'Missing required elements. '
                                               "Found {'C', 'O', 'N'}, need "
                                               "{'P', 'C', 'O', 'N'}"},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Missing required elements. '
                                               "Found {'C', 'O'}, need {'P', "
                                               "'C', 'O', 'N'}"},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Missing required elements. '
                                               "Found {'O', 'N', 'S', 'C', "
                                               "'Na'}, need {'P', 'C', 'O', "
                                               "'N'}"}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. The main issue seems to be that the code is failing to '
                 'properly identify valid structures because of the element '
                 "check. Looking at the error messages, it's saying molecules "
                 'are missing Hydrogen atoms, but this is likely because '
                 "RDKit's SMILES parser doesn't always include implicit "
                 'hydrogens in the atom list.\n'
                 '\n'
                 '2. The stereochemistry check might be too strict - we need '
                 'to account for different ways SMILES can represent the same '
                 'stereochemistry.\n'
                 '\n'
                 '3. We should also consider that some SMILES might represent '
                 'the molecule in ionic form (with [NH3+] and [O-] on the '
                 'phosphate).\n'
                 '\n'
                 "Here's an improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 27,
    'num_false_positives': 100,
    'num_true_negatives': 67276,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.2125984251968504,
    'recall': 1.0,
    'f1': 0.35064935064935066,
    'accuracy': 0.9985163865109862,
    'negative_predictive_value': 1.0}