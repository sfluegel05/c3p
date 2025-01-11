"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
"""
Classifies: glycosaminoglycan
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a glycosaminoglycan, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count key atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    
    if o_count < 3:  # Relaxed oxygen requirement
        return False, "Too few oxygen atoms for aminomonosaccharide structure"
    
    if n_count == 0:  # Still need nitrogen for amino groups
        return False, "No nitrogen atoms found - required for aminosugar"

    # Look for various sugar-like ring patterns
    patterns = [
        "[CR1]1[CR1][CR1][CR1][CR1]O1",  # pyranose
        "[CR1]1[CR1][CR1][CR1]O1",       # furanose
        "[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1][CR0,CR1]O1",  # modified pyranose
        "[CR0,CR1]1[CR0,CR1][CR0,CR1][CR0,CR1]O1"            # modified furanose
    ]
    
    total_sugar_matches = 0
    for pattern in patterns:
        sugar_pattern = Chem.MolFromSmarts(pattern)
        if sugar_pattern:
            matches = len(mol.GetSubstructMatches(sugar_pattern))
            total_sugar_matches += matches

    # Look for amino groups in various contexts
    amino_patterns = [
        "[NX3H,NX3H2][CH1,CH2][OH1,OR]",  # classic amino sugar
        "[NX3H,NX3H2][CR0,CR1][OR0,OR1]",  # modified amino sugar
        "[NX3H,NX3H2]C(=O)",               # N-acetyl group
        "[NX3]C[CR1]1O[CR1][CR1][CR1][CR1]1"  # N-substituted sugar
    ]
    
    total_amino_matches = 0
    for pattern in amino_patterns:
        amino_pattern = Chem.MolFromSmarts(pattern)
        if amino_pattern:
            matches = len(mol.GetSubstructMatches(amino_pattern))
            total_amino_matches += matches

    # Look for characteristic modifications
    sulfate = Chem.MolFromSmarts("OS(=O)(=O)[OH1,O-]")
    carboxyl = Chem.MolFromSmarts("C(=O)[OH1,O-]")
    acetyl = Chem.MolFromSmarts("NC(=O)C")
    
    has_sulfate = mol.HasSubstructMatch(sulfate) if sulfate else False
    has_carboxyl = mol.HasSubstructMatch(carboxyl) if carboxyl else False
    has_acetyl = mol.HasSubstructMatch(acetyl) if acetyl else False

    # Build features list
    features = []
    if total_sugar_matches > 0:
        features.append(f"Contains {total_sugar_matches} sugar-like rings")
    if total_amino_matches > 0:
        features.append(f"{total_amino_matches} amino groups")
    if has_sulfate:
        features.append("sulfate groups")
    if has_carboxyl:
        features.append("carboxyl groups")
    if has_acetyl:
        features.append("acetyl groups")

    # Classification criteria
    is_gag = (
        (total_sugar_matches > 0 or (o_count >= 4 and n_count >= 1)) and  # Either has sugar rings or enough O/N atoms
        total_amino_matches > 0 and                                        # Must have amino groups
        (has_sulfate or has_carboxyl or has_acetyl)                       # Must have at least one characteristic modification
    )

    if not features:
        return False, "No characteristic glycosaminoglycan features found"
    
    reason = ("Classified as glycosaminoglycan: " + ", ".join(features)) if is_gag else \
             ("Not classified as glycosaminoglycan: " + ", ".join(features))

    return is_gag, reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:18085',
                          'name': 'glycosaminoglycan',
                          'definition': 'Any polysaccharide containing a '
                                        'substantial proportion of '
                                        'aminomonosaccharide residues.',
                          'parents': ['CHEBI:22506'],
                          'xrefs': [   'KEGG:C02545',
                                       'Wikipedia:Glycosaminoglycan'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.038461538461538464 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](CO)O[C@H]1NC(=O)C[C@H](N)C(O)=O '
               'NAME: '
               'N(4)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)]-beta-D-Manp-(1->4)-beta-D-GlcpNAc-(1->4)-beta-D-GlcpNAc} '
               'REASON: CORRECT Classified as glycosaminoglycan: Contains 7 '
               'sugar rings, 1 amino sugar residues, 11 glycosidic linkages, '
               'carboxyl groups\n'
               ' * SMILES: '
               'CC(=O)N[C@@H]1[C@@H](O)[C@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H](CO)O[C@H]1NC(=O)C[C@H](N)C(O)=O '
               'NAME: '
               'N(4)-{beta-D-GlcNAc-(1->2)-alpha-D-Man-(1->3)-[beta-D-GlcNAc-(1->2)-alpha-D-Man-(1->6)]-beta-D-Man-(1->4)-beta-D-GlcNAc-(1->4)-beta-D-GlcNAc}-Asn '
               'REASON: CORRECT Classified as glycosaminoglycan: Contains 7 '
               'sugar rings, 1 amino sugar residues, 12 glycosidic linkages, '
               'carboxyl groups\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'ClC1=C2NC(=O)C=CC=CC=C[C@@H]([C@@H](O)CC(=O)C(C)=CC[C@@H](C=C[C@@H]([C@@H]([C@H](C=C(C(C=3C(C1=O)=C(C2=O)C=C(C)C3O)=O)C)C)O)C)O)C '
               'NAME: Naphthomycin B REASON: MISSED No sugar rings detected\n'
               ' * SMILES: '
               'S1S[C@@]23N4[C@@H]5[C@@H](OC(=O)C=6C=CC(=C(OC=7C=C([C@@H]([C@]1(N(C)C2=O)C4=O)O)C=CC7O)C6)OC)C=COC=C5[C@H]3O '
               'NAME: Emestrin REASON: MISSED No sugar rings detected\n'
               ' * SMILES: '
               'S1C2=NC(=C1)C(=O)N[C@@H](C=3SC=C(N3)C(=O)N[C@@H](C4=NC(C(N[C@@H]2C(C)C)=O)=C(O4)C)C)C '
               'NAME: Dendroamide A REASON: MISSED Molecular weight too low '
               'for glycosaminoglycan\n'
               ' * SMILES: '
               'S1C2=NC(=C1)C(=O)N[C@@H](C=3OC(C)=C(N3)C(=O)N[C@@H](C4=NC(C(N[C@@H]2CCSC)=O)=CS4)C)C '
               'NAME: Dendroamide B REASON: MISSED Too few oxygen atoms for '
               'polysaccharide structure\n'
               ' * SMILES: '
               'S([C@@]12N(C(=O)[C@@](SC)(CC=3C4=C(C=CC=C4)N(C3)[C@]56[C@H](N(COC1)C7=C5C=CC=C7)N8C(=O)[C@](SC)(N(C)C([C@@]8(SC)C6)=O)CO)N(C2=O)C)C)C '
               'NAME: Chaetocochin A REASON: MISSED No sugar rings detected\n'
               ' * SMILES: '
               'O=C1NCCC[C@@H]2NC(=O)C(C2=O)=C(O)C=CC=C[C@H]3[C@@H](C=CC=C1)C=C[C@H]4[C@H](C(=O)C)[C@@H](C[C@H]34)C '
               'NAME: Clifednamide F REASON: MISSED Molecular weight too low '
               'for glycosaminoglycan\n'
               ' * SMILES: '
               'O=C1C(=C[C@H](CC=C[C@@H]2[C@]3(C=4NC1=CC4)C(=O)N[C@H]([C@@H]3[C@H](C)C([C@H]2O)=C)CC=5C6=C(C=CC=C6)NC5)C)C '
               'NAME: Penochalasin C REASON: MISSED Too few oxygen atoms for '
               'polysaccharide structure\n'
               ' * SMILES: '
               'ClC1=C2NC(=O)C[C@H](O)C(=C[C@@H]([C@@H](O)CC(=O)[C@H](OC)C=CC=C(CC(=C2)C=C1O)C)C)C '
               'NAME: 19-Chloroproansamitocin REASON: MISSED Molecular weight '
               'too low for glycosaminoglycan\n'
               ' * SMILES: '
               'O1C2=C3C4=C5NC6(N=C5C(=NC(=O)C(=CC=C[C@H]([C@H](O)[C@@H]([C@H](O)[C@H]([C@@H](OC(=O)C)[C@@H]([C@H](O)C=CO[C@@]1(C3=O)C)C)C)C)C)C)C(O)=C4C(O)=C2C)CCN(CC6)CC(C)C '
               'NAME: 27-O-demethylrifabutin REASON: MISSED No sugar rings '
               'detected\n'
               ' * SMILES: '
               'O=C1O[C@H](C=CC(CCC(C=2C=3C(C=4[C@@]1(C=CC(=O)NC4C(=O)C3C=C(C)C2O)C)=O)=O)CC)C(O)C '
               'NAME: Hygrocin G REASON: MISSED No sugar rings detected\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O(C1=CC=2[C@]3([C@](N(CC3)C)(N(C2C=C1)C)[H])C)C(=O)N4CCC=5C(C4)=CC=CC5',
                                     'name': 'quilostigmine',
                                     'reason': 'Too few oxygen atoms for '
                                               'aminomonosaccharide structure'},
                                 {   'smiles': 'O[C@@H]1[C@]23[C@@]4(N(C[C@@]([C@]2(C[C@@]4([C@]56[C@]3(CC(=O)[C@](C5)(C([C@H]6O)=C)[H])[H])[H])[H])(CC1)C)CC)[H]',
                                     'name': 'Bullatine G',
                                     'reason': 'No characteristic '
                                               'glycosaminoglycan features '
                                               'found'},
                                 {   'smiles': 'NC1=NC=NC2=C1N=CN2[C@@H]3O[C@H](COP(=O)(O)O)[C@@H](OC(=O)[C@@H](N)CCC(O)=O)[C@H]3O',
                                     'name': "3'-L-glutamyl-AMP",
                                     'reason': 'Not classified as '
                                               'glycosaminoglycan: Contains 2 '
                                               'sugar-like rings, carboxyl '
                                               'groups'},
                                 {   'smiles': 'O1C2(C(C3(C(C4(C(CC3OC(=O)C)C(OC(=O)CC4)(C)C)C)CC2)C)CC15C6N(C=7C5=CC=CC7)C(=O)C(N6)C)C',
                                     'name': 'Teraspiridole C_130091',
                                     'reason': 'Not classified as '
                                               'glycosaminoglycan: acetyl '
                                               'groups'},
                                 {   'smiles': 'O(C1[C@@H](OC(=O)C)C(O[C@@H](OC2=C(OC3=C(C2=O)C(O)=CC(O[C@@H]4OC([C@@H](O)[C@H](O)C4O)CO)=C3CC=C(C)C)C5=CC=C(OC)C=C5)[C@H]1O)C)[C@@H]6OC[C@@H](O)[C@H](OC(=O)C)C6O',
                                     'name': 'Sempervirenoside A',
                                     'reason': 'No nitrogen atoms found - '
                                               'required for aminosugar'},
                                 {   'smiles': 'COC[C@]1(C(=O)C2CCN1CC2)CO',
                                     'name': '(2S)-2-(hydroxymethyl)-2-(methoxymethyl)-1-azabicyclo[2.2.2]octan-3-one',
                                     'reason': 'No characteristic '
                                               'glycosaminoglycan features '
                                               'found'},
                                 {   'smiles': 'Oc1c(C2CC(Cc3ccccc23)c2ccc(OCc3ccc(cc3)C(F)(F)F)cc2)c(=O)oc2ccccc12',
                                     'name': 'Flocoumafen',
                                     'reason': 'No nitrogen atoms found - '
                                               'required for aminosugar'},
                                 {   'smiles': 'O[C@H]1CC=2C(N(C=3C1=CC=CC3)C(=O)N)=CC=CC2',
                                     'name': '(S)-MHD',
                                     'reason': 'Too few oxygen atoms for '
                                               'aminomonosaccharide structure'},
                                 {   'smiles': 'S(OC=1C(O)=C(\\C=C\\C2=CC(O)=C(CC=C(C)C)C(O)=C2)C=CC1O)(O)(=O)=O',
                                     'name': '3-{(e)-2-[3,5-dihydroxy-4-(3-methyl-2-buten-1-yl)phenyl]vinyl}-2,6-dihydroxyphenyl '
                                             'hydrogen sulfate',
                                     'reason': 'No nitrogen atoms found - '
                                               'required for aminosugar'},
                                 {   'smiles': 'O1[C@@]2(O)[C@]3([C@](CC(=C[C@]3(C4=C1C=C5OC=6C(C5=C4C)=C(O)C=C(C6)CO)[H])C)(C7=C(OC=8C7=C(O)C=C(C8)CO)C2=O)C)[H]',
                                     'name': 'Asticolorin C',
                                     'reason': 'No nitrogen atoms found - '
                                               'required for aminosugar'}],
    'sample_false_negatives': [   {   'smiles': 'S1S[C@@]23N4[C@@H]5[C@@H](OC(=O)C=6C=CC(=C(OC=7C=C([C@@H]([C@]1(N(C)C2=O)C4=O)O)C=CC7O)C6)OC)C=COC=C5[C@H]3O',
                                      'name': 'Emestrin',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: acetyl '
                                                'groups'},
                                  {   'smiles': 'S1C2=NC(=C1)C(=O)N[C@@H](C=3SC=C(N3)C(=O)N[C@@H](C4=NC(C(N[C@@H]2C(C)C)=O)=C(O4)C)C)C',
                                      'name': 'Dendroamide A',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: 3 amino '
                                                'groups'},
                                  {   'smiles': 'S1C2=NC(=C1)C(=O)N[C@@H](C=3OC(C)=C(N3)C(=O)N[C@@H](C4=NC(C(N[C@@H]2CCSC)=O)=CS4)C)C',
                                      'name': 'Dendroamide B',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: 3 amino '
                                                'groups'},
                                  {   'smiles': 'S([C@@]12N(C(=O)[C@@](SC)(CC=3C4=C(C=CC=C4)N(C3)[C@]56[C@H](N(COC1)C7=C5C=CC=C7)N8C(=O)[C@](SC)(N(C)C([C@@]8(SC)C6)=O)CO)N(C2=O)C)C)C',
                                      'name': 'Chaetocochin A',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: acetyl '
                                                'groups'},
                                  {   'smiles': 'O=C1C(=C[C@H](CC=C[C@@H]2[C@]3(C=4NC1=CC4)C(=O)N[C@H]([C@@H]3[C@H](C)C([C@H]2O)=C)CC=5C6=C(C=CC=C6)NC5)C)C',
                                      'name': 'Penochalasin C',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: 1 amino '
                                                'groups, acetyl groups'},
                                  {   'smiles': 'O1C2=C3C4=C5NC6(N=C5C(=NC(=O)C(=CC=C[C@H]([C@H](O)[C@@H]([C@H](O)[C@H]([C@@H](OC(=O)C)[C@@H]([C@H](O)C=CO[C@@]1(C3=O)C)C)C)C)C)C)C(O)=C4C(O)=C2C)CCN(CC6)CC(C)C',
                                      'name': '27-O-demethylrifabutin',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: acetyl '
                                                'groups'},
                                  {   'smiles': 'S1C2=NC(=C1)C(=O)N[C@H](C=3OC(C)=C(N3)C(=O)NCC4=NC(C(N[C@@H]2CCSC)=O)=CS4)C',
                                      'name': 'Nostocyclamide M',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: 3 amino '
                                                'groups'},
                                  {   'smiles': 'S1[C@H](CC(=O)N[C@H](C=2SC[C@H](C(N[C@@H](C1)CO)=O)N2)CC3=CC=CC=C3)CCC(C)C',
                                      'name': 'Thioamycolamide D',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: 2 amino '
                                                'groups, acetyl groups'},
                                  {   'smiles': 'S1S[C@@]23N4[C@@H]5[C@@H](OC(=O)C=6C=CC(=C(OC=7C=C(C[C@]1(N(C)C2=O)C4=O)C=CC7O)C6)OC)C=COC=C5[C@H]3O',
                                      'name': 'Emestrin F',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: acetyl '
                                                'groups'},
                                  {   'smiles': 'S1SS[C@@]23N4[C@@H]5[C@@H](OC(=O)C=6C=CC(=C(OC=7C=C([C@@H](C(S1)(N(C)C2=O)C4=O)O)C=CC7OC)C6)OC)C=COC=C5[C@H]3O',
                                      'name': 'Emestrin E',
                                      'reason': 'Not classified as '
                                                'glycosaminoglycan: acetyl '
                                                'groups'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the previous outcomes and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Problem Analysis:\n'
                 '- The code correctly identified the complex polysaccharide '
                 'structures with multiple sugar rings and amino groups\n'
                 '- However, it failed to identify many valid '
                 'glycosaminoglycans that have different structural patterns\n'
                 '- The main issues appear to be:\n'
                 '  a) Too strict requirements for sugar rings (pyranose '
                 'pattern only)\n'
                 '  b) Too rigid amino sugar pattern matching\n'
                 '  c) Molecular weight and atom count thresholds may be too '
                 'restrictive\n'
                 '\n'
                 '2. Key Observations from False Negatives:\n'
                 '- Many compounds were missed due to "No sugar rings '
                 'detected" even though they contain sugar-like structures\n'
                 '- Some compounds were rejected due to molecular weight being '
                 'too low\n'
                 '- The pattern matching for amino sugars seems too specific\n'
                 '\n'
                 '3. Proposed Improvements:\n'
                 '- Broaden sugar ring detection to include more patterns '
                 '(furanose, modified rings)\n'
                 '- Relax the molecular weight threshold\n'
                 '- Add more flexible patterns for amino sugar detection\n'
                 '- Consider alternative indicators of glycosaminoglycan '
                 'structure\n'
                 '- Add detection of characteristic functional groups\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 203,
    'num_false_positives': 100,
    'num_true_negatives': 287,
    'num_false_negatives': 43,
    'num_negatives': None,
    'precision': 0.66996699669967,
    'recall': 0.8252032520325203,
    'f1': 0.7395264116575592,
    'accuracy': 0.7740916271721959,
    'negative_predictive_value': 0.8696969696969697}