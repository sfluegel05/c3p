"""
Classifies: CHEBI:33567 catecholamine
"""
"""
Classifies: catecholamine
Definition: 4-(2-Aminoethyl)pyrocatechol and derivatives formed by substitution
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_catecholamine(smiles: str):
    """
    Determines if a molecule is a catecholamine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a catecholamine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Convert to neutral form if possible (handle salt forms)
    mol = Chem.RemoveHs(mol)  # Remove explicit hydrogens
    
    # More flexible catechol patterns allowing different substitutions
    catechol_patterns = [
        Chem.MolFromSmarts("c1(O)c(O)cccc1"),  # 1,2-dihydroxy
        Chem.MolFromSmarts("c1c(O)c(O)ccc1"),   # 2,3-dihydroxy
        Chem.MolFromSmarts("c1cc(O)c(O)cc1"),   # 3,4-dihydroxy
        Chem.MolFromSmarts("c1c(O)cc(O)cc1"),   # 2,4-dihydroxy
        Chem.MolFromSmarts("c1(O)c(O)c(O)ccc1"), # Trihydroxy variants
        Chem.MolFromSmarts("c1(O)c(O)cc(O)cc1"),
    ]
    
    has_catechol = False
    for pat in catechol_patterns:
        if pat is not None and mol.HasSubstructMatch(pat):
            has_catechol = True
            break
            
    if not has_catechol:
        return False, "No catechol moiety found"
    
    # More comprehensive ethylamine patterns
    ethylamine_patterns = [
        # Basic patterns
        Chem.MolFromSmarts("[cR1]!@CC[NH2,NH1,NH0]"), # Direct ethylamine chain
        Chem.MolFromSmarts("[cR1]!@CC([*,H])N"),      # Allow substitution at beta carbon
        Chem.MolFromSmarts("[cR1]!@C([*,H])CN"),      # Allow substitution at alpha carbon
        # Hydroxylated variants
        Chem.MolFromSmarts("[cR1]!@CC(O)N"),          # Beta-hydroxyl
        Chem.MolFromSmarts("[cR1]!@C(O)CN"),          # Alpha-hydroxyl
        # Various amine substitutions
        Chem.MolFromSmarts("[cR1]!@CCN([*,H])[*,H]"), # Secondary amine
        Chem.MolFromSmarts("[cR1]!@CCN([*,H])([*,H])"), # Tertiary amine
    ]
    
    has_ethylamine = False
    for pat in ethylamine_patterns:
        if pat is not None and mol.HasSubstructMatch(pat):
            has_ethylamine = True
            break
            
    if not has_ethylamine:
        return False, "No ethylamine chain found"
    
    # Additional checks to avoid false positives
    
    # Check molecule size (catecholamines are relatively small)
    if mol.GetNumAtoms() > 40:  # Increased from 30 to allow for some larger derivatives
        return False, "Molecule too large for a typical catecholamine"
    
    # Check for maximum one benzene ring with catechol pattern
    aromatic_rings = mol.GetSubstructMatches(Chem.MolFromSmarts("a1aaaaa1"))
    if len(aromatic_rings) > 2:  # Allow max 2 rings for structures like dobutamine
        return False, "Too many aromatic rings for a catecholamine"
        
    # Check for presence of unlikely ring systems in catecholamines
    complex_ring_systems = [
        Chem.MolFromSmarts("C1=CC2=C(C=C1)C1=C(C=CC=C1)N2"), # Complex fused rings
        Chem.MolFromSmarts("C1CCCCC1"), # Cyclohexane
        Chem.MolFromSmarts("C1CCCC1"),  # Cyclopentane
    ]
    
    for ring in complex_ring_systems:
        if ring is not None and mol.HasSubstructMatch(ring):
            return False, "Contains ring systems not typical for catecholamines"
    
    return True, "Contains catechol moiety with ethylamine chain in appropriate arrangement"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:33567',
                          'name': 'catecholamine',
                          'definition': '4-(2-Aminoethyl)pyrocatechol '
                                        '[4-(2-aminoethyl)benzene-1,2-diol] '
                                        'and derivatives formed by '
                                        'substitution.',
                          'parents': [   'CHEBI:25375',
                                         'CHEBI:33566',
                                         'CHEBI:64365'],
                          'xrefs': ['KEGG:C02012'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.0392156862745098 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: CC1=C(C(=CC(=C1O)O)O)CCN NAME: '
               '5-(2-aminoethyl)-6-methylbenzene-1,2,4-triol REASON: CORRECT '
               'Contains catechol moiety with ethylamine chain\n'
               'False positives: SMILES: '
               'O1C=2C=3[C@@](NCCC3C(O)=C(OC)C2OC)(CC4=CC(OC5=CC=C(CC6=NCCC=7C6=CC1=C(OC)C7)C=C5)=C(OC)C=C4)[H] '
               'NAME: Thalsimidine REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               ' * SMILES: NC(Cc1cc(O)c(O)cc1O)C(O)=O NAME: 6-hydroxydopa '
               'REASON: WRONGLY CLASSIFIED Contains catechol moiety with '
               'ethylamine chain\n'
               ' * SMILES: COc1cc2C[C@@H]3N(C)CCc4ccc(OC)c(Oc2cc1OC)c34 NAME: '
               'Cularine REASON: WRONGLY CLASSIFIED Contains catechol moiety '
               'with ethylamine chain\n'
               ' * SMILES: '
               'O=C1N2C(C(=O)N[C@H](C(=O)O)CC(CO)CO)CC=3C(O)=C(O)C(=C(C3[C@H]2N[C@H]1CC=4C5=C(C=CC=C5)NC4)O)OC '
               'NAME: Imizoquin C REASON: WRONGLY CLASSIFIED Contains catechol '
               'moiety with ethylamine chain\n'
               ' * SMILES: '
               'O1[C@@H]2N([C@H]3C4=C(O)C(=C(OC)C(=C4[C@@H](CO)N5[C@H]3[C@@H]6[C@H]2C[C@H]([C@H]5C#N)N6C)O)C)CC1 '
               'NAME: Cyanocycline C REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               ' * SMILES: '
               'O=C1C(OC)=C(C(=O)C2=C1[C@@H](N3C[C@@H]4[C@@H](O)C=5C([C@H]([C@@H]3C2)N4C)=C(O)C(OC)=C(C)C5O)CNC(=O)C(=O)C)C '
               'NAME: Saframycin E REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               ' * SMILES: C[C@@H](C(=O)N[C@@H](C)C(=O)NCCC1=CC(=C(C=C1O)O)O)N '
               'NAME: '
               '(2S)-2-amino-N-[(2S)-1-oxo-1-[2-(2,4,5-trihydroxyphenyl)ethylamino]propan-2-yl]propanamide '
               'REASON: WRONGLY CLASSIFIED Contains catechol moiety with '
               'ethylamine chain\n'
               ' * SMILES: '
               'C=12C3=CC4=C(C1)OC=5C(O4)=C(C(=C6C5[C@@H](N(CC6)C)CC=7C=CC(=CC7)OC=8C=C(CC2=NCC3)C=CC8O)O)OC '
               'NAME: (+)-1,2-dehydrokohatine REASON: WRONGLY CLASSIFIED '
               'Contains catechol moiety with ethylamine chain\n'
               ' * SMILES: '
               'C[C@@H](C(=O)NCCC1=CC(=C(C=C1O)O)O)NC(=O)[C@H](C)NC(=O)OC(C)(C)C '
               'NAME: '
               'N-[(2S)-1-oxo-1-[[(2S)-1-oxo-1-[2-(2,4,5-trihydroxyphenyl)ethylamino]propan-2-yl]amino]propan-2-yl]carbamic '
               'acid tert-butyl ester REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               ' * SMILES: '
               'O(C1=C2C([N+](CCC2=CC=3OCOC13)(C)C)CC(=O)C4=C(OC)C5=C(OC=C5)C(OC)=C4O)C '
               'NAME: '
               '5-[2-(6-Hydroxy-4,7-dimethoxy-1-benzofuran-5-yl)-2-oxoethyl]-4-methoxy-6,6-dimethyl-5,6,7,8-tetrahydro[1,3]dioxolo[4,5-g]isoquinolin-6-ium '
               'REASON: WRONGLY CLASSIFIED Contains catechol moiety with '
               'ethylamine chain\n'
               ' * SMILES: '
               'O=C1C(OC)=C(C(=O)C2=C1[C@@H](N3[C@@H](O)[C@H]4[C@@H](OC)C=5C([C@@H]([C@@H]3C2)N4C)=C(O)C(OC)=C(C)C5O)CNC(=O)[C@@H](N)C)C '
               'NAME: Saframycin Mx1 REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               ' * SMILES: '
               'O=C1N([C@]23[C@H](C(C)(C)[C@@]4(C2)C(=O)NC=5C4=C(O)C=C6OC(C)(C)C=COC56)C[C@]17N(CC[C@@H]7C)C3)C '
               'NAME: SB203105 REASON: WRONGLY CLASSIFIED Contains catechol '
               'moiety with ethylamine chain\n'
               ' * SMILES: CN1CCc2ccc(O)c3Oc4cc5OCOc5cc4C[C@H]1c23 NAME: '
               'Cularicine REASON: WRONGLY CLASSIFIED Contains catechol moiety '
               'with ethylamine chain\n'
               ' * SMILES: '
               'COc1ccc2C[C@@H]3N(C)CCc4c(OC)c(OC)c(OC)c(Oc5cc6c(CCN=C6Cc6ccc(Oc1c2)cc6)cc5OC)c34 '
               'NAME: Thalsimine REASON: WRONGLY CLASSIFIED Contains catechol '
               'moiety with ethylamine chain\n'
               ' * SMILES: C1=C(C(=CC(=C1O)O)O)CC2C(=O)NC(=O)N2 NAME: '
               '5-[(2,4,5-trihydroxyphenyl)methyl]imidazolidine-2,4-dione '
               'REASON: WRONGLY CLASSIFIED Contains catechol moiety with '
               'ethylamine chain\n'
               ' * SMILES: '
               'O1C=2C=3C(N(CCC3C(OC)=C(OC)C2OC)C)CC4=CC(OC5=CC=C(CC6N(CCC=7C6=CC1=C(OC)C7)C)C=C5)=C(OC)C=C4 '
               'NAME: Isohernandezine REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               ' * SMILES: COc1cc2C[C@@H]3N(C)CCc4ccc(O)c(Oc2cc1OC)c34 NAME: '
               'Cularidine REASON: WRONGLY CLASSIFIED Contains catechol moiety '
               'with ethylamine chain\n'
               ' * SMILES: '
               'O=C1C(OC)=C(C(=O)C2=C1[C@@H](N3[C@@H](C#N)[C@@H]4CC=5C([C@H]([C@@H]3C2)N4C)=C(OC(=O)CO)C(OC)=C(C)C5O)CNC(=O)C(=O)C)C '
               'NAME: Saframycin R REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               ' * SMILES: '
               'COc1ccc2C[C@@H]3N(C)CCc4c(OC)c(OC)c(OC)c(Oc5cc6[C@H](Cc7ccc(Oc1c2)cc7)N(C)CCc6cc5OC)c34 '
               'NAME: Hernandezine REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               ' * SMILES: '
               'O=C1C(OC)=C(C(=O)C2=C1[C@@H](N3C[C@H]4C(=O)C=5C([C@@H]([C@@H]3C2)N4C)=C(O)C(OC)=C(C)C5O)CNC(=O)C(=O)C)C '
               'NAME: Saframycin D REASON: WRONGLY CLASSIFIED Contains '
               'catechol moiety with ethylamine chain\n'
               'False negatives: SMILES: '
               'C[C@@H](CCc1ccc(O)cc1)NCCc1ccc(O)c(O)c1 NAME: (S)-dobutamine '
               'REASON: MISSED No catechol (adjacent hydroxyl groups on '
               'benzene) found\n'
               ' * SMILES: CCCN(CCC)CCC1=CC(=C(C=C1)O)O NAME: '
               '4-[2-(dipropylamino)ethyl]benzene-1,2-diol REASON: MISSED No '
               'catechol (adjacent hydroxyl groups on benzene) found\n'
               ' * SMILES: O(C)C1=C(O)C=CC(C(O)CNC)=C1 NAME: Metanephrine '
               'REASON: MISSED No catechol (adjacent hydroxyl groups on '
               'benzene) found\n'
               ' * SMILES: C1=CC(=C(C(=C1CCN)O)O)O NAME: '
               '4-(2-aminoethyl)benzene-1,2,3-triol REASON: MISSED No catechol '
               '(adjacent hydroxyl groups on benzene) found\n'
               ' * SMILES: NC[C@H](O)c1ccc(O)c(O)c1 NAME: (R)-noradrenaline '
               'REASON: MISSED No catechol (adjacent hydroxyl groups on '
               'benzene) found\n'
               ' * SMILES: CNC[C@H](O)c1ccc(O)c(O)c1 NAME: (R)-adrenaline '
               'REASON: MISSED No catechol (adjacent hydroxyl groups on '
               'benzene) found\n'
               ' * SMILES: OC1=C(O)C=C(CCN)C=C1O NAME: 5-Hydroxydopamine '
               'REASON: MISSED No catechol (adjacent hydroxyl groups on '
               'benzene) found\n'
               ' * SMILES: NCC(O)c1ccc(O)c(O)c1 NAME: noradrenaline REASON: '
               'MISSED No catechol (adjacent hydroxyl groups on benzene) '
               'found\n'
               ' * SMILES: CC(N)C(O)c1ccc(O)c(O)c1 NAME: Nordephrine REASON: '
               'MISSED No catechol (adjacent hydroxyl groups on benzene) '
               'found\n'
               ' * SMILES: CC(C)NC[C@H](O)c1ccc(O)c(O)c1 NAME: L-isoprenaline '
               'REASON: MISSED No catechol (adjacent hydroxyl groups on '
               'benzene) found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No catechol moiety found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No catechol moiety found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No catechol moiety found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No catechol moiety found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No catechol moiety found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No catechol moiety found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No ethylamine chain found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No catechol moiety found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No catechol moiety found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No catechol moiety found'}],
    'sample_false_negatives': [],
    'attempt': 1,
    'reasoning': 'Let me analyze the issues with the previous code and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Pattern Matching Issues:\n'
                 '- The catechol patterns were too restrictive and missed '
                 'valid catechols with different substitution patterns\n'
                 "- The ethylamine patterns didn't capture all valid "
                 'variations of the aminoethyl chain\n'
                 '- The complete pattern check was too strict and rejected '
                 'valid structures\n'
                 '\n'
                 '2. False Positives:\n'
                 '- Many complex molecules containing fragments similar to '
                 'catecholamines were incorrectly classified\n'
                 '- Need better structural context checking\n'
                 '\n'
                 '3. False Negatives:\n'
                 '- Valid catecholamines like dopamine and noradrenaline were '
                 'missed\n'
                 '- Need to handle more substitution patterns and variations\n'
                 '\n'
                 "Here's an improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 30,
    'num_false_positives': 100,
    'num_true_negatives': 25692,
    'num_false_negatives': 0,
    'num_negatives': None,
    'precision': 0.23076923076923078,
    'recall': 1.0,
    'f1': 0.375,
    'accuracy': 0.9961273332816978,
    'negative_predictive_value': 1.0}