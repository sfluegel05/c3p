"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
"""
Classifies: CHEBI:24836 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for phosphate group(s)
    # More general phosphate pattern that matches mono-, di-, tri- phosphates
    phosphate_pattern = Chem.MolFromSmarts("[OX2]P(=O)([OX2,OX1])[OX2,OX1]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"
    
    # Look for 5'-phosphate connection - CH2-O-P
    sugar_phosphate = Chem.MolFromSmarts("[CH2X4]OP(~[OX1,OX2])~[OX1,OX2]")
    if not mol.HasSubstructMatch(sugar_phosphate):
        return False, "No 5'-phosphate group found"
        
    # Check for furanose ring with more flexible patterns
    sugar_ring_patterns = [
        # Basic furanose ring without stereochemistry
        Chem.MolFromSmarts("[CH2][C]1O[C][C][C]1"),
        # Ribose-like pattern (more general)
        Chem.MolFromSmarts("[CH2]([OH,O])[C]1O[C]([#6,#7])[C][C]1[OH,O]"),
        # Deoxyribose-like pattern (more general)
        Chem.MolFromSmarts("[CH2]([OH,O])[C]1O[C]([#6,#7])[C][C]1"),
        # Alternative sugar pattern
        Chem.MolFromSmarts("[CH2]OP~[OX1,OX2][C]1O[C][C][C]1")
    ]
    
    found_sugar = False
    for pattern in sugar_ring_patterns:
        if mol.HasSubstructMatch(pattern):
            found_sugar = True
            break
    
    if not found_sugar:
        return False, "No furanose ring found"

    # Nucleobase patterns (more general)
    purine_patterns = [
        # Basic purine scaffold (more flexible)
        Chem.MolFromSmarts("[#7,#6]1[#7,#6][#7,#6][#6]2[#7,#6]1[#7,#6][#7,#6][#7,#6]2"),
        # Modified purine
        Chem.MolFromSmarts("[#7,#6]1[#7,#6][#7,#6][#6]2[#7,#6]1[#7,#6][#6]([#7,O,S])[#7,#6]2")
    ]
    
    pyrimidine_patterns = [
        # Basic pyrimidine (more flexible)
        Chem.MolFromSmarts("[#7,#6]1[#7,#6][#7,#6][#6]([#7,O,S])[#7,#6][#7,#6]1"),
        # Alternative pyrimidine
        Chem.MolFromSmarts("[#7]1[#6][#6][#7][#6]([O,S,N])[#7]1")
    ]
    
    has_purine = any(mol.HasSubstructMatch(pattern) for pattern in purine_patterns)
    has_pyrimidine = any(mol.HasSubstructMatch(pattern) for pattern in pyrimidine_patterns)
    
    if not (has_purine or has_pyrimidine):
        return False, "No nucleobase (purine/pyrimidine) found"

    # Count phosphate groups
    p_o_p_pattern = Chem.MolFromSmarts("[PX4]~[OX2]~[PX4]")
    p_o_p_count = len(mol.GetSubstructMatches(p_o_p_pattern))
    phosphate_count = len(mol.GetSubstructMatches(phosphate_pattern))
    
    if p_o_p_count > 3 or phosphate_count > 4:
        return False, "Too many phosphate groups"
        
    # Determine phosphate type
    if p_o_p_count == 0:
        phos_type = "mono"
    elif p_o_p_count == 1:
        phos_type = "di"
    elif p_o_p_count == 2:
        phos_type = "tri"
    else:
        phos_type = "tetra"
    
    base_type = "purine" if has_purine else "pyrimidine"
    return True, f"Found {phos_type}-phosphorylated nucleoside with {base_type} base"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16701',
                          'name': "nucleoside 5'-phosphate",
                          'definition': 'A ribosyl or deoxyribosyl derivative '
                                        'of a pyrimidine or purine base in '
                                        'which C-5 of the ribose ring is '
                                        'mono-, di-, tri- or '
                                        'tetra-phosphorylated.',
                          'parents': ['CHEBI:29075'],
                          'xrefs': ['KEGG:C01117'],
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
               'C[n+]1cn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c2nc(N)[nH]c(=O)c12 '
               "NAME: 7-methylguanosine 5'-phosphate(1+) REASON: MISSED No "
               'ribose/deoxyribose ring found\n'
               ' * SMILES: '
               'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CC4=CNC5=C4C=CC=C5)[C@H]3O '
               "NAME: 3'-L-tryptophyl-AMP REASON: MISSED No ribose/deoxyribose "
               'ring found\n'
               ' * SMILES: '
               'Nc1ncnc2n([C@H]3C[C@H](O)[C@@H](COP(O)(=O)OP(O)(O)=O)O3)c(O)nc12 '
               'NAME: 8-hydroxy-dADP REASON: MISSED No ribose/deoxyribose ring '
               'found\n'
               ' * SMILES: Nc1ccn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)n1 '
               "NAME: 2'-deoxycytosine 5'-monophosphate REASON: MISSED No "
               'ribose/deoxyribose ring found\n'
               ' * SMILES: '
               'Cc1nc2c(N)ncnc2n1[C@@H]1O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]1O '
               "NAME: 8-methyladenosine 5'-monophosphate REASON: MISSED No "
               'ribose/deoxyribose ring found\n'
               ' * SMILES: '
               'O1[C@@H]([C@H]([C@H]([C@@H]1N2C(NC(C=C2)=O)=O)O)OC(=C)C(=O)O)COP(=O)(O)O '
               "NAME: 3'-(enolpyruvyl)uridine 5'-monophosphate REASON: MISSED "
               'No ribose/deoxyribose ring found\n'
               ' * SMILES: '
               '[C@H]1([C@H](O[C@H]([C@@H]1O)N2C(NC3=C2N=C(NC3=O)N)=O)COP(OP(O)(=O)O)(=O)O)O '
               'NAME: 8-oxo-GDP REASON: MISSED No ribose/deoxyribose ring '
               'found\n'
               ' * SMILES: '
               'Cc1cc2N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(O)=O)C3=NC(=O)NC(=O)C3N3CCC(C)(C)c(c1C)c23 '
               'NAME: prenyl-FMNH2 REASON: MISSED No ribose/deoxyribose ring '
               'found\n'
               ' * SMILES: '
               'Nc1ncnc2ncn([C@@H]3O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]3O)c12 '
               'NAME: 7-(5-phospho-alpha-D-ribosyl)adenine REASON: MISSED No '
               'ribose/deoxyribose ring found\n'
               ' * SMILES: '
               'CC(=O)Nc1ccn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c(=O)n1 '
               "NAME: N(4)-acetylcytidine 5'-monophosphate REASON: MISSED No "
               'ribose/deoxyribose ring found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'S(=O)(=O)(C1=C2C(=CC3=C1NC=4C=CC=CC34)[C@@]5([C@H]([C@](C(=O)O)([C@@H](O)CC5)C)CC2)C)C6=CC7=C(NC8=C7C=C9[C@@]%10([C@H]([C@](C(=O)O)([C@@H](O)CC%10)C)CCC9=C8)C)C=C6',
                                     'name': 'Sulfadixiamycin C',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'CNC(O)=O',
                                     'name': 'methylcarbamic acid',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'CCNC(=O)NC1=CC2=C(C=C1)OC[C@H]3[C@@H](CC[C@H](O3)CC(=O)N[C@@H](C)C4=CC=CC=C4)N(C2=O)C',
                                     'name': '2-[(2S,4aR,12aR)-8-(ethylcarbamoylamino)-5-methyl-6-oxo-2,3,4,4a,12,12a-hexahydropyrano[2,3-c][1,5]benzoxazocin-2-yl]-N-[(1S)-1-phenylethyl]acetamide',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'O([C@H]1O[C@@H]([C@@H](O[C@@H]2O[C@@H]([C@@H](O[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO[C@]4(O[C@H]([C@H](NC(=O)C)[C@@H](O)C4)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]2NC(=O)C)CO)[C@H](O)[C@@H]1O[C@@H]5O[C@@H]([C@@H](O[C@@H]6O[C@@H]([C@H](O)[C@H](O)[C@H]6O)CO[C@]7(O[C@H]([C@H](NC(=O)C)[C@@H](O)C7)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]5NC(=O)C)CO)CO)[C@H]8[C@H](O)[C@H](O[C@@H](O[C@H]9[C@H](O)[C@@H](NC(=O)C)[C@@H](O[C@@H]9CO)O[C@H]%10[C@H](O)[C@@H](NC(=O)C)C(O[C@@H]%10CO[C@@H]%11O[C@H]([C@@H](O)[C@@H](O)[C@@H]%11O)C)O)[C@H]8O)CO[C@H]%12O[C@@H]([C@@H](O)[C@H](O)[C@@H]%12O[C@@H]%13O[C@@H]([C@@H](O[C@@H]%14O[C@@H]([C@H](O)[C@H](O[C@]%15(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%15)[C@H](O)[C@H](O)CO)C(O)=O)[C@H]%14O)CO)[C@H](O)[C@H]%13NC(=O)C)CO)CO[C@@H]%16O[C@@H]([C@@H](O[C@@H]%17O[C@@H]([C@H](O)[C@H](O)[C@H]%17O)CO[C@]%18(O[C@H]([C@H](NC(=O)C)[C@@H](O)C%18)[C@H](O)[C@H](O)CO)C(O)=O)[C@H](O)[C@H]%16NC(=O)C)CO',
                                     'name': 'CID 91851985',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'O(C(=O)C(C1C(CN2C(C1)C=3NC=4C(C3CC2)=CC=CC4)CC)=COC)C',
                                     'name': 'Methyl '
                                             '2-(3-ethyl-1,2,3,4,6,7,12,12b-octahydroindolo[2,3-a]quinolizin-2-yl)-3-methoxyprop-2-enoate',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'O[C@H](/C=C/C=C/C=C/[C@H](O)[C@H](O)C=C)[C@H](O)/C=C/C',
                                     'name': 'Separacene C',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'C[C@@H]1O[C@@H](O[C@@H]2[C@@H](CO)O[C@@H](O[C@@H]3[C@@H](O)C(O)O[C@H](CO)[C@@H]3O)[C@H](NC(C)=O)[C@H]2O[C@@H]2O[C@H](CO)[C@H](O)[C@H](OS(O)(=O)=O)[C@H]2O)[C@@H](O)[C@H](O)[C@@H]1O',
                                     'name': 'alpha-L-Fucp-(1->4)-[beta-D-Galp3S-(1->3)]-beta-D-GlcpNAc-(1->3)-D-Galp',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'C1=CC=CC2=C1C(N([C@H](C(N2)=O)CC=3C=CC(=CC3)OC)C)=O',
                                     'name': "(S)-4'-methoxycyclopeptine",
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'O=C(N[C@@H](CC=1C=2C(NC1)=CC=CC2)C(O)=O)[C@@H](NC(=O)[C@@H](N)C(C)C)C(C)C',
                                     'name': 'Val-Val-Trp',
                                     'reason': 'No phosphate group found'},
                                 {   'smiles': 'C=1C(=C(C=CC1/C=C/CO)OC(CO)C(O)C=2C=C(C(=CC2)O)OC)OC',
                                     'name': 'guaiacylglycerol beta-coniferyl '
                                             'ether',
                                     'reason': 'No phosphate group found'}],
    'sample_false_negatives': [   {   'smiles': 'O1[C@@H]([C@H]([C@H]([C@@H]1N2C(NC(C=C2)=O)=O)O)OC(=C)C(=O)O)COP(=O)(O)O',
                                      'name': "3'-(enolpyruvyl)uridine "
                                              "5'-monophosphate",
                                      'reason': 'No nucleobase '
                                                '(purine/pyrimidine) found'},
                                  {   'smiles': 'Cc1cc2N(C[C@H](O)[C@H](O)[C@H](O)COP(O)(O)=O)C3=NC(=O)NC(=O)C3N3CCC(C)(C)c(c1C)c23',
                                      'name': 'prenyl-FMNH2',
                                      'reason': 'No furanose ring found'},
                                  {   'smiles': 'Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)([O-])=O)O2)c(=O)[nH]c1=O',
                                      'name': 'dTMP(-)',
                                      'reason': 'No nucleobase '
                                                '(purine/pyrimidine) found'},
                                  {   'smiles': 'NCc1cn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c(=S)[nH]c1=O',
                                      'name': '5-aminomethyl-2-thiouridine '
                                              "5'-monophosphate",
                                      'reason': 'No nucleobase '
                                                '(purine/pyrimidine) found'},
                                  {   'smiles': 'Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(=O)OP(O)(O)=O)O2)c(=O)[nH]c1=O',
                                      'name': 'dTDP',
                                      'reason': 'No nucleobase '
                                                '(purine/pyrimidine) found'},
                                  {   'smiles': 'C[C@@H](OP(O)(=O)OC[C@@H](O)[C@@H](O)[C@@H](O)Cn1c2cc(C)c(C)cc2nc2c1nc(=O)[nH]c2=O)[C@H](N)C(O)=O',
                                      'name': 'FMN-L-threonine',
                                      'reason': 'No furanose ring found'},
                                  {   'smiles': 'O[C@H]1[C@@H](O)[C@@H](O[C@@H]1COP(O)(O)=O)n1cc(F)c(=O)[nH]c1=O',
                                      'name': '5-fluorouridine '
                                              "5'-monophosphate",
                                      'reason': 'No nucleobase '
                                                '(purine/pyrimidine) found'},
                                  {   'smiles': 'C1(C(C(=NC(N1[C@@H]2O[C@H](COP(O)(=O)O)[C@H](C2)O)=O)N)O)O',
                                      'name': "2'-deoxy-5,6-dihydroxy-5,6-dihydrocytidine "
                                              "5'-monophosphate",
                                      'reason': 'No nucleobase '
                                                '(purine/pyrimidine) found'},
                                  {   'smiles': 'Cn1c(=N)ccn([C@@H]2O[C@H](COP(O)(O)=O)[C@@H](O)[C@H]2O)c1=O',
                                      'name': 'N(3)-methylcytidine '
                                              "5'-monophosphate",
                                      'reason': 'No nucleobase '
                                                '(purine/pyrimidine) found'},
                                  {   'smiles': 'O[C@H]1C[C@@H](O[C@@H]1COP(O)(=O)OP(O)(O)=O)n1ccc(=O)[nH]c1=O',
                                      'name': 'dUDP',
                                      'reason': 'No nucleobase '
                                                '(purine/pyrimidine) found'}],
    'attempt': 3,
    'reasoning': "Let's analyze what went wrong and how to improve:\n"
                 '\n'
                 '1. The main issue seems to be with the sugar ring detection. '
                 'The current SMARTS patterns for ribose/deoxyribose are too '
                 'specific and rigid, causing false negatives.\n'
                 '\n'
                 '2. Looking at the failed examples:\n'
                 "   - They all contain valid nucleoside 5'-phosphate "
                 'structures\n'
                 '   - The sugar ring patterns are failing to match despite '
                 'correct structures\n'
                 '   - The stereochemistry specifications in the SMARTS '
                 'patterns might be too restrictive\n'
                 '\n'
                 '3. Improvements needed:\n'
                 '   - Simplify the sugar ring pattern to be more flexible\n'
                 '   - Make stereochemistry optional in the patterns\n'
                 '   - Focus on the essential connectivity rather than exact '
                 'configuration\n'
                 '   - Add more variations of sugar-base connections\n'
                 '   - Consider using recursive SMARTS for more complex '
                 'patterns\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 123,
    'num_false_positives': 100,
    'num_true_negatives': 6999,
    'num_false_negatives': 43,
    'num_negatives': None,
    'precision': 0.5515695067264574,
    'recall': 0.7409638554216867,
    'f1': 0.6323907455012855,
    'accuracy': 0.9803165863730213,
    'negative_predictive_value': 0.9938937801760863}