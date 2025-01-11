"""
Classifies: CHEBI:31998 N-acylphytosphingosine
"""
"""
Classifies: N-acylphytosphingosine
Definition: A ceramide that is phytosphingosine having a fatty acyl group attached to the nitrogen.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_N_acylphytosphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylphytosphingosine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        tuple: (bool, str) - (is_N_acylphytosphingosine, reason)
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for amide group (N-C(=O))
    amide_pattern = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No amide group found"

    # Look for phytosphingosine core with more flexible pattern
    # Allow for variations in connectivity while maintaining core structure
    phyto_core = Chem.MolFromSmarts("[OX2H1,OX2H0][#6]-[#6]([NX3])-[#6]([OX2H1,OX2H0])-[#6]([OX2H1,OX2H0])-[#6]")
    if not mol.HasSubstructMatch(phyto_core):
        return False, "No phytosphingosine core structure found"

    # Count hydroxyl groups (both free and substituted)
    oh_pattern = Chem.MolFromSmarts("[OX2H1,OX2H0]-[#6]")
    oh_matches = len(mol.GetSubstructMatches(oh_pattern))
    if oh_matches < 3:
        return False, f"Insufficient hydroxyl groups (found {oh_matches}, need ≥3)"

    # Check for fatty acyl chain - more flexible pattern
    # Look for carbon chain attached to amide
    fatty_chain = Chem.MolFromSmarts("[NX3][CX3](=[OX1])[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(fatty_chain):
        return False, "No suitable fatty acyl chain found"

    # Basic size checks
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18:
        return False, f"Insufficient carbon atoms (found {c_count}, need ≥18)"

    mol_weight = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_weight < 300:
        return False, f"Molecular weight too low ({mol_weight:.1f} Da)"

    # Check for possible sugar modifications
    sugar_pattern = Chem.MolFromSmarts("[#6]1-[#8]-[#6]-[#6]-[#6]-[#6]1")
    has_sugar = mol.HasSubstructMatch(sugar_pattern)

    # Additional checks for chain length
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain length too short for N-acylphytosphingosine"

    # Check nitrogen count
    n_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_count < 1:
        return False, "No nitrogen found"
    if n_count > 3 and not has_sugar:
        return False, f"Too many nitrogens for basic structure (found {n_count})"

    # Check for long chain characteristic of sphingolipids
    long_chain = Chem.MolFromSmarts("[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "Missing characteristic long carbon chain"

    return True, "Contains phytosphingosine backbone with N-acyl group and required hydroxyl groups"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:31998',
                          'name': 'N-acylphytosphingosine',
                          'definition': 'A ceramide that is phytosphingosine '
                                        'having a fatty acyl group attached to '
                                        'the nitrogen.',
                          'parents': ['CHEBI:139051'],
                          'xrefs': [   'KEGG:C12145',
                                       'LIPID_MAPS_instance:LMSP01030000'],
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
               '[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)OC[C@@H]([C@@H]([C@@H](CCCCCCCCCCCCCC)O)O)NC(CCCCCCCCCCC=2C=CC=CC2)=O '
               'NAME: '
               '1-O-(alpha-D-galactopyranosyl)-N-(11-phenylundecanoyl)phytosphingosine '
               'REASON: MISSED No phytosphingosine core structure found\n'
               ' * SMILES: '
               'O([C@@H]1[C@@H]([C@H]([C@H]([C@H](O1)CO)OCC=2C=CC(=CC2)C)O)O)C[C@@H]([C@@H]([C@@H](CCCCCCCCCCCCCC)O)O)NC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O '
               'NAME: '
               '1-O-[4-O-(4-methylbenzyl)-alpha-D-galactosyl]-N-hexacosanoylphytosphingosine '
               'REASON: MISSED No phytosphingosine core structure found\n'
               ' * SMILES: '
               'O=C(N[C@H]([C@H](O)[C@H](O)CCCCCCCCCCCCCC)CO)[C@H](O)CCCCCCCCCCCCCCCCCCCCC '
               "NAME: (2S, 3S, 4R)-2-[(2'R)-2'-hydroxytricosanoylamino]-1, 3, "
               '4-octadecanetriol REASON: MISSED No phytosphingosine core '
               'structure found\n'
               ' * SMILES: '
               'C12(CCO[C@@H]3[C@@H]([C@H]([C@H](O[C@@H]3CO)OC[C@@H]([C@@H]([C@@H](CCCCCCCCCCCCCC)O)O)NC(CCCCCCCCCCCCCCCCCCCCCCCCC)=O)O)O)C[C@@H]4CC(C1)C[C@H](C4)C2 '
               'NAME: '
               '1-O-[4-O-(2-adamantylethyl)-alpha-D-galactosyl]-N-hexacosanoylphytosphingosine '
               'REASON: MISSED No phytosphingosine core structure found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCC[C@@H](O)[C@@H](O)[C@H](CO[C@H]1O[C@H](CO)[C@H](O)[C@H](O)[C@H]1O)NC(=O)CCCCCCCC '
               'NAME: '
               '1-O-(alpha-D-galactopyranosyl)-N-nonanoylphytosphingosine '
               'REASON: MISSED No phytosphingosine core structure found\n'
               ' * SMILES: '
               '[C@H]1([C@@H]([C@H]([C@H]([C@H](O1)CO)O)O)O)OC[C@@H]([C@@H]([C@@H](CCCCCCCCCCCCCC)O)O)NC(CCCCCCCC2=CC=C(C=C2)C=3C=CC=CC3)=O '
               'NAME: '
               "N-[8-([1,1'-biphenyl]-4-yl)octanoyl]-1-O-(alpha-D-galactopyranosyl)phytosphingosine "
               'REASON: MISSED No phytosphingosine core structure found\n'
               ' * SMILES: '
               'OC(C(O)C(NC(=O)C(O)CCCCCCCCCCCCCCCCCCCCCCCC)CO)CCCCCCCCCCCCCC '
               'NAME: '
               'N-(2R-Hydroxyhexacosanoyl)-2S-amino-1,3S,4R-octadecanetriol '
               'REASON: MISSED No phytosphingosine core structure found\n'
               ' * SMILES: '
               'OC(C(O)C(NC(=O)C(O)CCCCCCCCCCCCCCCCCCCC)CO)CCCCCCCCCCCCCC '
               'NAME: 2-Hydroxy-N-(1,3,4-trihydroxy-2-octadecanyl)docosanamide '
               'REASON: MISSED No phytosphingosine core structure found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC '
               'NAME: N-triacontanoyphytosphingosine REASON: MISSED No '
               'phytosphingosine core structure found\n'
               ' * SMILES: '
               'CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)[C@H](O)CCCCCCCCCCCCCC '
               'NAME: N-octadecanoyl-4-hydroxysphinganine REASON: MISSED No '
               'phytosphingosine core structure found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'O=C1C(=C2C=C3[C@]([C@@H](C(C)C)[C@@H]([C@H]3O)OC(=O)C)(C)CC[C@]2(C)CC1)COC(=O)C',
                                     'name': 'Dahliane E',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'No phytosphingosine core '
                                               'structure found'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'No phytosphingosine core '
                                               'structure found'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'No phytosphingosine core '
                                               'structure found'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'No amide group found'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'No amide group found'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'No suitable fatty acyl chain '
                                               'found'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'No amide group found'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'No phytosphingosine core '
                                               'structure found'}],
    'sample_false_negatives': [   {   'smiles': 'O[C@@H](CCCCCC/C=C/[C@H](O)[C@@H](O)[C@H](O)[C@H](NC(=O)C)C(O)=O)CCCCCC',
                                      'name': 'Sphingofungin D',
                                      'reason': 'No suitable fatty acyl chain '
                                                'found'}],
    'attempt': 3,
    'reasoning': 'Let me analyze what went wrong and propose improvements:\n'
                 '\n'
                 '1. The main issue seems to be that the SMARTS pattern for '
                 'the phytosphingosine core was too restrictive:\n'
                 '   - The previous pattern '
                 '`[CH2][OX2]-[CH]([NX3])-[CH]([OX2])-[CH]([OX2])-[CH2]` was '
                 'requiring specific connectivity and hybridization states\n'
                 '   - Many valid structures were missed because they have '
                 'slight variations in connectivity or atom types\n'
                 '\n'
                 '2. The hydroxyl group checks were too specific:\n'
                 '   - We should focus on the overall pattern rather than '
                 'counting specific types of OH groups\n'
                 '   - The structure should have 3 hydroxyl groups and one '
                 "amine group that's acylated\n"
                 '\n'
                 '3. The fatty acyl chain check was too restrictive:\n'
                 '   - The previous pattern required exactly 6 carbons in a '
                 'specific arrangement\n'
                 '   - We should allow for more flexibility in chain length '
                 'and structure\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 93,
    'num_false_positives': 100,
    'num_true_negatives': 38089,
    'num_false_negatives': 1,
    'num_negatives': None,
    'precision': 0.48186528497409326,
    'recall': 0.9893617021276596,
    'f1': 0.6480836236933797,
    'accuracy': 0.9973617532586265,
    'negative_predictive_value': 0.9999737463901286}