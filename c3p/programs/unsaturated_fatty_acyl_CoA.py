"""
Classifies: CHEBI:51006 unsaturated fatty acyl-CoA
"""
"""
Classifies: unsaturated fatty acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_unsaturated_fatty_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acyl-CoA based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an unsaturated fatty acyl-CoA, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for adenine using multiple possible patterns
    adenine_patterns = [
        "c1ncnc2[nH]cnc12",  # Basic adenine core
        "c1ncnc2ncnc12",     # Alternative representation
        "c1nc(N)nc2[nH]cnc12",  # With amino group
        "c1nc(N)nc2ncnc12",     # Another amino form
        "[nH]1cnc2c(ncnc2n1)",  # Different tautomer
        "n1cnc2c(N)ncnc12"      # Yet another form
    ]
    
    found_adenine = False
    for pattern in adenine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_adenine = True
            break
    
    if not found_adenine:
        return False, "No adenine moiety found"

    # Look for phosphate groups characteristic of CoA
    phosphate_patterns = [
        "OP(O)(=O)OP(O)(=O)O",
        "P(O)(O)(=O)OP(O)(O)=O"
    ]
    
    found_phosphates = False
    for pattern in phosphate_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_phosphates = True
            break
            
    if not found_phosphates:
        return False, "Missing characteristic phosphate groups of CoA"
    
    # Check for pantetheine arm with thioester
    pantetheine_patterns = [
        "NCCC(=O)NCCS",
        "SCCNC(=O)CCNC"
    ]
    
    found_pantetheine = False
    for pattern in pantetheine_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_pantetheine = True
            break
            
    if not found_pantetheine:
        return False, "Missing characteristic pantetheine arm of CoA"

    # Check for thioester linkage
    thioester_patterns = [
        "C(=O)S",
        "SC(=O)"
    ]
    
    found_thioester = False
    for pattern in thioester_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_thioester = True
            break
            
    if not found_thioester:
        return False, "No thioester linkage found"

    # Count number of double bonds in the fatty acid portion
    # Look for C=C bonds that are not part of the adenine ring
    double_bond_pattern = Chem.MolFromSmarts("C=C")
    double_bonds = mol.GetSubstructMatches(double_bond_pattern)
    
    # Get aromatic atoms
    aromatic_atoms = {atom.GetIdx() for atom in mol.GetAromaticAtoms()}
    
    # Filter out double bonds where either carbon is aromatic
    non_aromatic_double_bonds = [
        bond for bond in double_bonds 
        if not (bond[0] in aromatic_atoms or bond[1] in aromatic_atoms)
    ]
    
    if len(non_aromatic_double_bonds) == 0:
        return False, "No carbon-carbon double bonds found in fatty acid portion"

    # Check for reasonable chain length (at least 4 carbons in fatty acid portion)
    carbon_chain_patterns = [
        "CCCC",
        "C~C~C~C"  # More flexible pattern allowing any bonds between carbons
    ]
    
    found_chain = False
    for pattern in carbon_chain_patterns:
        if mol.HasSubstructMatch(Chem.MolFromSmarts(pattern)):
            found_chain = True
            break
            
    if not found_chain:
        return False, "Fatty acid portion too short"

    # Success case - molecule has all required features
    return True, f"Contains CoA thioester linkage and {len(non_aromatic_double_bonds)} C=C double bonds in fatty acid portion"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51006',
                          'name': 'unsaturated fatty acyl-CoA',
                          'definition': 'A fatty acyl-CoA that results from '
                                        'the formal condensation of the thiol '
                                        'group of coenzyme A with the carboxy '
                                        'group of any unsaturated fatty acid.',
                          'parents': ['CHEBI:37554'],
                          'xrefs': ['PMID:13152086'],
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
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (6Z,9Z,12Z,15Z,18Z,21Z)-3-oxotetracosahexaenoyl-CoA '
               'REASON: MISSED No adenine moiety found\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: '
               '(3R,18Z,21Z,24Z,27Z,30Z,33Z)-3-hydroxyhexatriacontahexaenoyl-CoA '
               'REASON: MISSED No adenine moiety found\n'
               ' * SMILES: '
               'CCCC\\C=C/C\\C=C/C\\C=C/CCCCCC(=O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (9Z,12Z,15Z)-3-oxoicosatrienoyl-CoA REASON: MISSED No '
               'adenine moiety found\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (3R,13Z,16Z,19Z,22Z)-3-hydroxyoctacosatetraenoyl-CoA '
               'REASON: MISSED No adenine moiety found\n'
               ' * SMILES: '
               'CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCCC\\C=C\\C(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (2E,19Z,22Z,25Z,28Z)-tetratriacontapentaenoyl-CoA '
               'REASON: MISSED No adenine moiety found\n'
               ' * SMILES: '
               'CC(=C)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: 3-methylbut-3-enoyl-CoA REASON: MISSED No adenine moiety '
               'found\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCCCCCCCC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (15Z,18Z,21Z,24Z,27Z)-triacontapentaenoyl-CoA REASON: '
               'MISSED No adenine moiety found\n'
               ' * SMILES: '
               'CCCCCCCC\\C=C/CCCCCCCCCCC[C@@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: (3R,15Z)-3-hydroxytetracosenoyl-CoA REASON: MISSED No '
               'adenine moiety found\n'
               ' * SMILES: '
               'CC(=CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12)C(C)(C)C=CC(O)=O '
               'NAME: 3,4,4-trimethylhepta-2,5-dienoyl-CoA REASON: MISSED No '
               'adenine moiety found\n'
               ' * SMILES: '
               'CC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CC[C@H](O)CC(=O)SCCNC(=O)CCNC(=O)[C@H](O)C(C)(C)COP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1OP(O)(O)=O)n1cnc2c(N)ncnc12 '
               'NAME: '
               '(3S,6Z,9Z,12Z,15Z,18Z,21Z)-3-hydroxytetracosahexaenoyl-CoA '
               'REASON: MISSED No adenine moiety found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No adenine moiety found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No adenine moiety found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No adenine moiety found'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No adenine moiety found'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No adenine moiety found'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No adenine moiety found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No adenine moiety found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No adenine moiety found'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'Missing characteristic '
                                               'phosphate groups of CoA'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No adenine moiety found'}],
    'sample_false_negatives': [],
    'attempt': 2,
    'reasoning': "Let's analyze what went wrong and how to improve the code:\n"
                 '\n'
                 '1. The main issue seems to be that the adenine pattern '
                 'matching is failing for valid structures. Looking at the '
                 'error messages, all molecules were rejected due to "No '
                 'adenine moiety found" despite clearly containing the adenine '
                 'group.\n'
                 '\n'
                 '2. The problem likely lies in the SMARTS patterns being too '
                 'specific or rigid. We need to make them more flexible to '
                 'account for different representations and tautomeric forms '
                 'of the adenine moiety.\n'
                 '\n'
                 '3. The CoA structure can be represented in different ways in '
                 'SMILES notation, and we need to ensure our patterns can '
                 'match all valid forms.\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 258,
    'num_false_positives': 100,
    'num_true_negatives': 29514,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.7206703910614525,
    'recall': 1.0,
    'f1': 0.8376623376623377,
    'accuracy': 0.9966523835029459,
    'negative_predictive_value': 1.0}