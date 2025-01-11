"""
Classifies: CHEBI:51129 nitrohydrocarbon
"""
"""
Classifies: CHEBI:51753 nitrohydrocarbon
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_nitrohydrocarbon(smiles: str):
    """
    Determines if a molecule is a nitrohydrocarbon based on its SMILES string.
    A nitrohydrocarbon is a hydrocarbon where one or more hydrogens are replaced by nitro groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nitrohydrocarbon, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of nitro groups (-NO2)
    # More specific SMARTS pattern for nitro groups
    nitro_pattern = Chem.MolFromSmarts('[NX3+](=[OX1])([O-])-[#6]')
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    if not nitro_matches:
        return False, "No nitro groups found"

    # Get all atoms
    atoms = mol.GetAtoms()
    
    # Check if molecule contains only C, H, N, and O
    allowed_atoms = {1, 6, 7, 8}  # H, C, N, O atomic numbers
    for atom in atoms:
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains elements other than C, H, N, O"

    # Check for carbon atoms (must have at least one)
    carbon_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 6]
    if not carbon_atoms:
        return False, "No carbon atoms found"

    # Check for unwanted functional groups that would make it not a pure nitrohydrocarbon
    unwanted_groups = [
        ('[CX3](=O)[OX2H1]', 'carboxylic acid'),  # More specific SMARTS
        ('[CX3](=O)[O-]', 'carboxylate'),
        ('[CX3]=O', 'ketone/aldehyde'),
        ('[NX3;H2,H1;!$(NC=O)]', 'amine'),  # More specific, excludes nitro
        ('[SX2]', 'thiol/sulfide'),
        ('[Cl]', 'chloride'),
        ('[Br]', 'bromide'),
        ('[F]', 'fluoride'),
        ('[I]', 'iodide'),
        ('[OX2H]', 'hydroxyl'),  # More specific SMARTS
        ('[OX2H]-[CX4]', 'alcohol'),  # More specific SMARTS
    ]
    
    for pattern, group_name in unwanted_groups:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol and mol.HasSubstructMatch(pattern_mol):
            return False, f"Contains {group_name} group"

    # Count nitro groups
    num_nitro = len(nitro_matches)
    
    # Verify all nitrogen atoms are part of nitro groups
    n_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 7]
    if len(n_atoms) != len(nitro_matches):
        return False, "Contains nitrogen atoms not in nitro groups"
    
    # Verify all oxygen atoms are part of nitro groups
    o_atoms = [atom for atom in atoms if atom.GetAtomicNum() == 8]
    if len(o_atoms) != 2 * len(nitro_matches):
        return False, "Contains oxygen atoms not in nitro groups"
    
    return True, f"Hydrocarbon with {num_nitro} nitro group{'s' if num_nitro > 1 else ''}"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:51129',
                          'name': 'nitrohydrocarbon',
                          'definition': 'A C-nitro compound that is a '
                                        'hydrocarbon in which one or more of '
                                        'the hydrogens has been replaced by '
                                        'nitro groups.',
                          'parents': ['CHEBI:35716'],
                          'xrefs': ['PMCID:PMC3909368'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: SMILES: '
               'OC=1C=2C=3N(CCCC3C=CC2)C(=O)C1[N+]([O-])=O NAME: '
               '7-hydroxy-6-nitro-2,3-dihydro-1H,5H-pyrido[3,2,1-ij]quinolin-5-one '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: '
               'Cc1c(c(C)c(c(c1[N+]([O-])=O)C(C)(C)C)[N+]([O-])=O)[N+]([O-])=O '
               'NAME: musk xylene REASON: WRONGLY CLASSIFIED Hydrocarbon with '
               '3 nitro groups\n'
               ' * SMILES: [O-][N+](=O)c1ccc(CCc2ccc(cc2)[N+]([O-])=O)cc1 '
               'NAME: 1,2-Bis(4-nitrophenyl)ethane REASON: WRONGLY CLASSIFIED '
               'Hydrocarbon with 2 nitro groups\n'
               ' * SMILES: CC[C@H](C)c1cc(cc(c1O)[N+]([O-])=O)[N+]([O-])=O '
               'NAME: (S)-dinoseb REASON: WRONGLY CLASSIFIED Hydrocarbon with '
               '2 nitro groups\n'
               ' * SMILES: C(=C/C1=CC=C(N=N1)N)\\C=2OC(=CC2)[N+]([O-])=O NAME: '
               'nifurprazine REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 '
               'nitro group\n'
               ' * SMILES: C1=C(C=C(C(=C1[N+](=O)[O-])O)O)C=C(C#N)C#N NAME: '
               '2-[(3,4-dihydroxy-5-nitrophenyl)methylidene]propanedinitrile '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: '
               'CC1=CN=NN1C2=CC=C(C=C2)C(=NNC3=CC=C(C=C3)[N+](=O)[O-])C NAME: '
               'N-[1-[4-(5-methyl-1-triazolyl)phenyl]ethylideneamino]-4-nitroaniline '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: C1CC(=NNC2=CC=CC=C2[N+](=O)[O-])C3=CC=CC=C31 NAME: '
               'N-(2,3-dihydroinden-1-ylideneamino)-2-nitroaniline REASON: '
               'WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: '
               'O=C1N(C(=O)N=C2N(CCCCCCCC)C3=C(C=C21)C=C([N+]([O-])=O)C=C3)C4=CC=CC=C4 '
               'NAME: CAY10625 REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 '
               'nitro group\n'
               ' * SMILES: '
               'C1=CC(=CC=C1/C=C(\\C)/C=C(\\C)/C=C(\\C)/C=C(\\C)/CCC2=C(C(C(=C(O2)[O-])C)=O)C)[N+]([O-])=O '
               'NAME: demethyldeoxyspectinabilin(1-) REASON: WRONGLY '
               'CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: '
               'C1C=CC2C1CC2=NNC3=C(C=C(C=C3)[N+](=O)[O-])[N+](=O)[O-] NAME: '
               'N-(6-bicyclo[3.2.0]hept-3-enylideneamino)-2,4-dinitroaniline '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 2 nitro groups\n'
               ' * SMILES: '
               'O1C(=NN=C1C2=CC([N+]([O-])=O)=CC=C2)C(C3=CC=C(CC(C)C)C=C3)C '
               'NAME: '
               '2-[1-(4-isobutylphenyl)ethyl]-5-(3-nitrophenyl)-1,3,4-oxadiazole '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: C1=CC(=C(C=C1C=C(C#N)C#N)O)[N+](=O)[O-] NAME: '
               '2-[(3-hydroxy-4-nitrophenyl)methylidene]propanedinitrile '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: [O-][N+](=O)C=1N(C(=NC1)C(C)C)C NAME: Ipronidazole '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: [O-][N+](=O)c1ccc(cc1)\\C=C\\c1ccccc1 NAME: '
               '(E)-4-nitrostilbene REASON: WRONGLY CLASSIFIED Hydrocarbon '
               'with 1 nitro group\n'
               ' * SMILES: '
               'C1=CC(=CC=C1/C=C(\\C)/C=C(\\C)/C=C(\\C)/C=C(\\C)/CCC2=C(C(C(=C(O2)O)C)=O)C)[N+]([O-])=O '
               'NAME: demethyldeoxyspectinabilin REASON: WRONGLY CLASSIFIED '
               'Hydrocarbon with 1 nitro group\n'
               ' * SMILES: '
               'C1CCC(=NNC2=C(C=C(C=C2)[N+](=O)[O-])[N+](=O)[O-])C(C1)C3=CC=CC=N3 '
               'NAME: '
               '2,4-dinitro-N-[[2-(2-pyridinyl)cyclohexylidene]amino]aniline '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 2 nitro groups\n'
               ' * SMILES: '
               'CC(=NNC1=C(C=C(C=C1)[N+](=O)[O-])[N+](=O)[O-])CC(C2=CC=CC=C2)C3=C(C4=CC=CC=C4OC3=O)O '
               'NAME: '
               '3-[3-[(2,4-dinitrophenyl)hydrazinylidene]-1-phenylbutyl]-4-hydroxy-1-benzopyran-2-one '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 2 nitro groups\n'
               ' * SMILES: '
               'C1=CC(=CC=C1/C=C(/C=C(/CCC2=C(C(C(=C(O2)[O-])C)=O)C)\\C)\\C)[N+]([O-])=O '
               'NAME: demethylluteothin(1-) REASON: WRONGLY CLASSIFIED '
               'Hydrocarbon with 1 nitro group\n'
               ' * SMILES: CC(C)(C)c1cc(cc(c1O)[N+]([O-])=O)[N+]([O-])=O NAME: '
               'Dinoterb REASON: WRONGLY CLASSIFIED Hydrocarbon with 2 nitro '
               'groups\n'
               ' * SMILES: '
               'C1=CC2=C(C(=C1)[N+](=O)[O-])N=C(C=C2)C=CC3=CC=C(C=C3)[N+](=O)[O-] '
               'NAME: 8-nitro-2-[2-(4-nitrophenyl)ethenyl]quinoline REASON: '
               'WRONGLY CLASSIFIED Hydrocarbon with 2 nitro groups\n'
               ' * SMILES: OC1=C(C=C([N+]([O-])=O)C=C1)\\C=N\\C=2C=C(C=CC2)C#C '
               'NAME: 2-{[(3-eth-1-ynylphenyl)imino]methyl}-4-nitrophenol '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: [O-][N+](=O)c1ccc(cc1)\\C=C\\c1ccc(cc1)[N+]([O-])=O '
               "NAME: 4,4'-Dinitrostilbene REASON: WRONGLY CLASSIFIED "
               'Hydrocarbon with 2 nitro groups\n'
               ' * SMILES: '
               'CC(C)(C)C1=CC=C(C=C1)CN2C=C(N=C2)C3=CC=C(C=C3)[N+](=O)[O-] '
               'NAME: '
               '1-[(4-tert-butylphenyl)methyl]-4-(4-nitrophenyl)imidazole '
               'REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 nitro group\n'
               ' * SMILES: [O-][N+](=O)c1ccc(\\C=C(/c2ccco2)c2n[nH]c(=N)o2)o1 '
               'NAME: Furamizole REASON: WRONGLY CLASSIFIED Hydrocarbon with 1 '
               'nitro group\n'
               'False negatives: SMILES: '
               '[O-][N+](=O)c1ccc2ccc(cc2c1)[N+]([O-])=O NAME: '
               '2,7-dinitronaphthalene REASON: MISSED No connected carbon '
               'framework found\n'
               ' * SMILES: [O-][N+](=O)c1ccc([N+]([O-])=O)c2ccccc12 NAME: '
               '1,4-dinitronaphthalene REASON: MISSED No connected carbon '
               'framework found\n'
               ' * SMILES: [O-][N+](=O)c1ccccc1 NAME: nitrobenzene REASON: '
               'MISSED No connected carbon framework found\n'
               ' * SMILES: [O-][N+](=O)c1ccc2ccccc2c1 NAME: 2-nitronaphthalene '
               'REASON: MISSED No connected carbon framework found\n'
               ' * SMILES: '
               'C1=CC2=C(C=CC3=C2C(=C1)C(=O)N(C3=O)CCC(=O)O)[N+](=O)[O-] NAME: '
               '3-(6-nitro-1,3-dioxo-2-benzo[de]isoquinolinyl)propanoic acid '
               'REASON: MISSED Contains carboxylic acid group\n'
               ' * SMILES: CC(C)[N+]([O-])=O NAME: 2-nitropropane REASON: '
               'MISSED Contains amine group\n'
               ' * SMILES: CCC[N+]([O-])=O NAME: 1-nitropropane REASON: MISSED '
               'Contains amine group\n'
               ' * SMILES: [O-][N+](=O)c1cccc2cccc([N+]([O-])=O)c12 NAME: '
               '1,8-dinitronaphthalene REASON: MISSED No connected carbon '
               'framework found\n'
               ' * SMILES: [O-][N+](=O)C1=C(C(=C(C(=C1[2H])[2H])[2H])[2H])[2H] '
               'NAME: nitrobenzene-d5 REASON: MISSED No connected carbon '
               'framework found\n'
               ' * SMILES: [O-][N+](=O)C=C NAME: nitroethene REASON: MISSED '
               'Contains amine group\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'C1C=C(N=C2N1NN=N2)C3=CC=CC=C3',
                                     'name': '5-phenyl-1,7-dihydrotetrazolo[1,5-a]pyrimidine',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': 'O(C(=O)CCCCCCCCCCC/C=C\\CCCCCCCC)[C@H](COC(=O)CCCCCCCCC/C=C\\CCCCCC)COC(=O)CCCCCCC/C=C\\CCCC',
                                     'name': 'TG(14:1(9Z)/22:1(13Z)/18:1(11Z))',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': 'CC(=O)N[C@H]1[C@H](O)O[C@H](CO)[C@@H](O[C@@H]2O[C@H](CO)[C@@H](O[C@@H]3O[C@H](CO[C@H]4O[C@H](CO[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H](O)[C@H](O[C@H]5O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]5O)[C@@H]4O)[C@@H](O)[C@H](O[C@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@@H]4O[C@@H]4O[C@H](CO)[C@@H](O)[C@H](O)[C@H]4NC(C)=O)[C@@H]3O)[C@H](O)[C@H]2NC(C)=O)[C@@H]1O',
                                     'name': 'beta-D-GlcpNAc-(1->2)-alpha-D-Manp-(1->3)-{alpha-D-Manp-(1->3)-[alpha-D-Manp-(1->6)]-alpha-D-Manp-(1->6)}-beta-D-Manp-(1->4)-beta-GlcpNAc-(1->4)-beta-D-GlcpNAc',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': '[C@@]12([C@@]([C@]([C@@H](CC1)C)(CCC(CC)C)C)(CCC[C@@H]2C)[H])C',
                                     'name': 'clerodane',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': 'S(O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@@H]([C@@H](O)[C@H](O)CO[C@]2(O[C@H]([C@H](NC(=O)C)[C@@H](O)C2)[C@H](O)[C@H](O)CO)C(O)=O)[C@@H](NC(=O)C)CO)O[C@@H]1CO)(O)(=O)=O',
                                     'name': '(2R,4S,5R,6R)-5-Acetamido-2-[(2R,3S,4R,5S)-5-acetamido-4-[(2R,3R,4R,5R,6R)-3,4-dihydroxy-6-(hydroxymethyl)-5-sulfooxyoxan-2-yl]oxy-2,3,6-trihydroxyhexoxy]-4-hydroxy-6-[(1R,2R)-1,2,3-trihydroxypropyl]oxane-2-carboxylic '
                                             'acid',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': 'O=C(NC(CC(C)C)C(=O)N[C@@H](CCCN=C(N)N)C=O)C(NC(=O)CC)CC(C)C',
                                     'name': 'Leupeptin Pr-LL',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': 'COC1=CC=C(C=C1)CCNC(=O)C(=C2C3=CC=CC=C3C(=N2)NC(=O)C4=CC=CS4)C#N',
                                     'name': 'N-[3-[1-cyano-2-[2-(4-methoxyphenyl)ethylamino]-2-oxoethylidene]-1-isoindolyl]-2-thiophenecarboxamide',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': 'O=C1O[C@@H](CC=2C1=C(OC)C(O)=CC2)CCC[C@@H](O)C',
                                     'name': 'Penicimarin C',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': 'CO[C@@H]1[C@@H]2[C@@H](C[C@@H]3[C@@](O2)([C@H]([C@@H]([C@H](O3)C(=O)OC)O)O)O)OC1N4C=NC5=C4N=CN=C5N',
                                     'name': 'LSM-4497',
                                     'reason': 'No nitro groups found'},
                                 {   'smiles': 'C[C@@H]1[C@H](O)CCC2=CC[C@H](C[C@]12C)C(C)=C',
                                     'name': '1-deoxycapsidiol',
                                     'reason': 'No nitro groups found'}],
    'sample_false_negatives': [   {   'smiles': 'C1=CC2=C(C=CC3=C2C(=C1)C(=O)N(C3=O)CCC(=O)O)[N+](=O)[O-]',
                                      'name': '3-(6-nitro-1,3-dioxo-2-benzo[de]isoquinolinyl)propanoic '
                                              'acid',
                                      'reason': 'Contains carboxylic acid '
                                                'group'},
                                  {   'smiles': 'O\\N=C1\\C(=O)Nc2c1cc(c1CCCCc21)N(=O)=O',
                                      'name': 'NS-102',
                                      'reason': 'Contains ketone/aldehyde '
                                                'group'},
                                  {   'smiles': 'C1=CC=C2C(=C1)C=CC(=C2O)NN=C3C4=C(C=C(C=C4)[N+](=O)[O-])C(=CC3=O)S(=O)(=O)O',
                                      'name': '4-[(1-hydroxy-2-naphthalenyl)hydrazinylidene]-7-nitro-3-oxo-1-naphthalenesulfonic '
                                              'acid',
                                      'reason': 'Contains elements other than '
                                                'C, H, N, O'},
                                  {   'smiles': 'C1=CC(=CC(=C1)N2C(=O)C3=CC(=CC4=CC(=CC(=C43)C2=O)[N+](=O)[O-])[N+](=O)[O-])C(=O)O',
                                      'name': '3-(5,8-dinitro-1,3-dioxo-2-benzo[de]isoquinolinyl)benzoic '
                                              'acid',
                                      'reason': 'Contains carboxylic acid '
                                                'group'},
                                  {   'smiles': '[O-][N+](=O)c1cccc2C3OC3C=Cc12',
                                      'name': '1-Nitronaphthalene-5,6-oxide',
                                      'reason': 'Contains oxygen atoms not in '
                                                'nitro groups'},
                                  {   'smiles': '[O-][N+](=O)c1cccc2C=CC3OC3c12',
                                      'name': '1-Nitronaphthalene-7,8-oxide',
                                      'reason': 'Contains oxygen atoms not in '
                                                'nitro groups'}],
    'attempt': 1,
    'reasoning': 'Let me analyze the issues with the previous code:\n'
                 '\n'
                 '1. The main problem seems to be with the SMARTS pattern for '
                 'nitro groups and carbon framework detection. Many valid '
                 'nitrohydrocarbons are being missed.\n'
                 '\n'
                 '2. The carbon framework detection is too strict - it '
                 'requires C~C which misses single-carbon compounds.\n'
                 '\n'
                 '3. Some of the unwanted group checks are too broad and may '
                 'be catching legitimate parts of nitrohydrocarbons.\n'
                 '\n'
                 "4. The amine check with 'CN' is incorrectly flagging nitro "
                 'groups as amines.\n'
                 '\n'
                 "Here's an improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 20,
    'num_false_positives': 34,
    'num_true_negatives': 142240,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.37037037037037035,
    'recall': 0.7692307692307693,
    'f1': 0.5,
    'accuracy': 0.9997189037245257,
    'negative_predictive_value': 0.9999578195520437}