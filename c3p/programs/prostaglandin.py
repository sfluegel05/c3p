"""
Classifies: CHEBI:26333 prostaglandin
"""
"""
Classifies: prostaglandin compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prostaglandin(smiles: str):
    """
    Determines if a molecule is a prostaglandin based on its SMILES string.
    Prostaglandins are naturally occurring C20 compounds derived from prostanoic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prostaglandin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbons - should be approximately 20 (allow variation for derivatives)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 15 or c_count > 35:
        return False, f"Carbon count ({c_count}) outside typical range for prostaglandins"

    # More flexible cyclopentane core pattern that includes various oxidation states
    # and possible peroxide bridges (as in PGH series)
    core_patterns = [
        "[#6]1~[#6]~[#6]~[#6]~[#6]1",  # Basic 5-membered ring
        "[#6]1~[#6]~[#6](~[#8,#6])~[#6]~[#6]1",  # With substituents
        "[#6]1~[#6]2~[#6]~[#8]~[#8]~[#6]2~[#6]~[#6]~[#6]1"  # Peroxide bridge pattern
    ]
    
    ring_found = False
    for pattern in core_patterns:
        core_pat = Chem.MolFromSmarts(pattern)
        if core_pat and mol.HasSubstructMatch(core_pat):
            ring_found = True
            break
    
    if not ring_found:
        return False, "Missing characteristic cyclopentane ring structure"

    # Check for characteristic chains with double bonds
    # More flexible pattern that captures both cis and trans configurations
    chain_pattern = Chem.MolFromSmarts("[#6]~[#6]=[#6]~[#6]")
    if not mol.HasSubstructMatch(chain_pattern):
        return False, "Missing characteristic alkenyl chains"

    # Check for carboxylic acid or derivatives (including esters and amides)
    acid_patterns = [
        "[CX3](=[OX1])[OX2H]",  # Carboxylic acid
        "[CX3](=[OX1])[OX2][#6]",  # Ester
        "[CX3](=[OX1])[NX3]",  # Amide
        "[CX3](=[OX1])[O-]"  # Carboxylate
    ]
    
    acid_found = False
    for pattern in acid_patterns:
        acid_pat = Chem.MolFromSmarts(pattern)
        if acid_pat and mol.HasSubstructMatch(acid_pat):
            acid_found = True
            break
            
    if not acid_found:
        return False, "Missing carboxylic acid or derivative group"

    # Count oxygen atoms (allow wider range for derivatives)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, f"Insufficient oxygen atoms ({o_count}) for prostaglandin"

    # Look for oxygen substituents (more flexible pattern)
    oxygen_patterns = [
        "[#6]~[#8H1]",  # Hydroxyl
        "[#6]=[#8]",    # Ketone
        "[#8]~[#8]"     # Peroxide
    ]
    
    oxygen_found = False
    for pattern in oxygen_patterns:
        o_pat = Chem.MolFromSmarts(pattern)
        if o_pat and mol.HasSubstructMatch(o_pat):
            oxygen_found = True
            break
            
    if not oxygen_found:
        return False, "Missing characteristic oxygen substituents"

    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 600:  # Widened range for derivatives
        return False, f"Molecular weight ({mol_wt:.1f}) outside typical range"

    return True, "Matches prostaglandin structural features: cyclopentane ring, oxygen substituents, characteristic chains"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26333',
                          'name': 'prostaglandin',
                          'definition': 'Naturally occurring compounds derived '
                                        'from the parent C20 acid, prostanoic '
                                        'acid.',
                          'parents': ['CHEBI:26347'],
                          'xrefs': [   'CAS:11000-26-3',
                                       'LIPID_MAPS_class:LMFA0301'],
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
               'C1[C@H]2[C@@H]([C@H]([C@@H]1OO2)/C=C/[C@H](CCCC(C)O)O)CCCCCCC(O)=O '
               'NAME: 19-hydroxyprostaglandin H1 REASON: MISSED Error in '
               'oxygen SMARTS pattern\n'
               ' * SMILES: '
               'C(O)CCCCCC[C@@H]1[C@H]([C@@H](C[C@@H]1O)O)/C=C/[C@H](CCCCC)O '
               'NAME: prostaglandin F1alpha alcohol REASON: MISSED No '
               'carboxylic acid or derivative found\n'
               ' * SMILES: '
               'CCCCC[C@H](O)\\C=C\\C1=CCC(=O)[C@@H]1C\\C=C/CCCC(O)=O NAME: '
               'prostaglandin C2 REASON: MISSED No cyclopentane ring found\n'
               ' * SMILES: '
               'CCCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)C[C@H](O)[C@@H]1CC(=O)CCCCC(O)=O '
               'NAME: 6-oxoprostaglandin F1alpha REASON: MISSED Error in '
               'oxygen SMARTS pattern\n'
               ' * SMILES: '
               'CCCCCC(=O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1C\\C=C/CCCC(O)=O '
               'NAME: 15-dehydro-prostaglandin E2 REASON: MISSED No '
               'cyclopentane ring found\n'
               ' * SMILES: '
               'CCCCCCCC(=O)CC[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CCCC(O)=O '
               'NAME: unoprostone REASON: MISSED Error in oxygen SMARTS '
               'pattern\n'
               ' * SMILES: '
               'CC(C)OC(=O)CCC\\C=C/C[C@H]1[C@@H](O)C[C@@H](O)[C@@H]1\\C=C\\C(F)(F)COc1ccccc1 '
               'NAME: tafluprost REASON: MISSED No carboxylic acid or '
               'derivative found\n'
               ' * SMILES: '
               'C(\\C[C@H]1[C@H](CC([C@@H]1/C=C/[C@H](CCCCC)O)=O)O)=C\\CCCC(=O)OC(CO)CO '
               'NAME: prostaglandin D2 2-glyceryl ester REASON: MISSED No '
               'cyclopentane ring found\n'
               ' * SMILES: '
               'CCCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)CC(=O)[C@@H]1C\\C=C/CCCC(O)=O '
               'NAME: prostaglandin E2 REASON: MISSED No cyclopentane ring '
               'found\n'
               ' * SMILES: '
               'C(\\[C@H]1[C@@H](CC([C@@H]1C/C=C\\CCCC(O)=O)=O)O)=C/[C@H](CCCCCO)O '
               'NAME: 20-hydroxyprostaglandin E2 REASON: MISSED No '
               'cyclopentane ring found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'Cc1cc(C)c2Cc3c(C)cc(C)c(Cc4c(C)cc(C)c(Cc5c(C)cc(C)c(Cc1c2C)c5C)c4C)c3C',
                                     'name': '4,6,10,12,16,18,22,24,25,26,27,28-dodecamethylcalix[4]arene',
                                     'reason': 'Carbon count (40) outside '
                                               'typical range for '
                                               'prostaglandins'},
                                 {   'smiles': 'O[C@@H]([C@H](NC(=O)[C@@H](N)CCC(O)=O)C(=O)N[C@@H](CC(C)C)C(O)=O)C',
                                     'name': 'Glu-Thr-Leu',
                                     'reason': 'Missing characteristic '
                                               'cyclopentane ring structure'},
                                 {   'smiles': 'CCOc1ccc(NC(=O)C(C)O)cc1',
                                     'name': 'p-Lactophenetide',
                                     'reason': 'Carbon count (11) outside '
                                               'typical range for '
                                               'prostaglandins'},
                                 {   'smiles': 'O=C1NCC=2C1=C3C(N([C@@H]4O[C@H]([C@@H](O)[C@H]([C@H]4OC)O)C)C5=C3C=CC=C5)=C6NC7=C(C26)C=CC=C7',
                                     'name': "3'-epi-5'-methoxy-K252d",
                                     'reason': 'Missing characteristic '
                                               'cyclopentane ring structure'},
                                 {   'smiles': 'O=C(OC1C(O)C(OC(C1O)CO)OCC2OC(OCC(OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\C/C=C\\C/C=C\\CC)C(O)C(C2O)O)CCCCCCCCCCCCCCC',
                                     'name': '[(2S)-2-hexadecanoyloxy-3-[(2S,3R,4S,5S,6R)-6-[[(2S,3R,4S,5S,6R)-4-hexadecanoyloxy-3,5-dihydroxy-6-(hydroxymethyl)tetrahydropyran-2-yl]oxymethyl]-3,4,5-trihydroxy-tetrahydropyran-2-yl]oxy-propyl] '
                                             '(9E,12E,15E)-octadeca-9,12,15-trienoate',
                                     'reason': 'Carbon count (65) outside '
                                               'typical range for '
                                               'prostaglandins'},
                                 {   'smiles': 'O=C1C2=C(O)C=C(OC)C=C2C(=O)C3=C1[C@@H]([C@@H](O)[C@]([C@@H]3O)(O)C)C',
                                     'name': 'Altersolanol G',
                                     'reason': 'Missing characteristic '
                                               'cyclopentane ring structure'},
                                 {   'smiles': '[H][C@]1(O[C@](O)(C[C@H](O)[C@H]1NC(=O)CO)C([O-])=O)[C@H](O)[C@H](O)CO',
                                     'name': 'N-glycoloyl-alpha-neuraminate',
                                     'reason': 'Carbon count (11) outside '
                                               'typical range for '
                                               'prostaglandins'},
                                 {   'smiles': 'OC(C(O)C/C=C\\C/C=C\\C/C=C\\CC)C/C=C\\C/C=C\\CCC(O)=O',
                                     'name': '10,11-DiHDPE',
                                     'reason': 'Missing characteristic '
                                               'cyclopentane ring structure'},
                                 {   'smiles': '[Na+].[H][C@]12SCC(C)=C(N1C(=O)[C@H]2NC(=O)[C@H](N)c1ccccc1)C([O-])=O',
                                     'name': 'cephalexin sodium',
                                     'reason': 'Missing characteristic '
                                               'cyclopentane ring structure'},
                                 {   'smiles': 'C1=CC=C(C=C1)COC2=CC=CC=C2C=NN3C(=NN=N3)N',
                                     'name': '1-[(2-phenylmethoxyphenyl)methylideneamino]-5-tetrazolamine',
                                     'reason': 'Missing characteristic '
                                               'cyclopentane ring structure'}],
    'sample_false_negatives': [   {   'smiles': 'C(O)CCCCCC[C@@H]1[C@H]([C@@H](C[C@@H]1O)O)/C=C/[C@H](CCCCC)O',
                                      'name': 'prostaglandin F1alpha alcohol',
                                      'reason': 'Missing carboxylic acid or '
                                                'derivative group'},
                                  {   'smiles': 'CCCCC[C@H](O)\\C=C\\[C@H]1[C@H](O)C[C@H](O)[C@@H]1C\\C=C/CCCCN(C)C',
                                      'name': 'prostaglandin F2alpha '
                                              'dimethylamine',
                                      'reason': 'Missing carboxylic acid or '
                                                'derivative group'},
                                  {   'smiles': 'CCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O',
                                      'name': 'nonacosanoic acid',
                                      'reason': 'Missing characteristic '
                                                'cyclopentane ring structure'}],
    'attempt': 3,
    'reasoning': 'Let me analyze the issues with the previous code and propose '
                 'improvements:\n'
                 '\n'
                 '1. Core Structure Recognition:\n'
                 '- The previous code was too strict about the cyclopentane '
                 'ring pattern. Prostaglandins can have various substituents '
                 'and oxidation states on the ring.\n'
                 '- Some prostaglandins (like PGH series) have peroxide '
                 'bridges in the ring.\n'
                 '- Need to make the ring pattern more flexible to account for '
                 'these variations.\n'
                 '\n'
                 '2. Chain Recognition:\n'
                 '- The chain pattern "CC=CC" was too simple and missed many '
                 'valid structures.\n'
                 '- Need to better capture the characteristic α and ω chains '
                 'with their various modifications.\n'
                 '- Should account for both cis and trans double bonds.\n'
                 '\n'
                 '3. Oxygen Substituents:\n'
                 '- The oxygen pattern was too specific and caused false '
                 'negatives.\n'
                 '- Need to handle various oxidation states (hydroxy, keto, '
                 'peroxy) more flexibly.\n'
                 '\n'
                 '4. Carboxylic Acid Recognition:\n'
                 '- The acid pattern was good but should also consider various '
                 'derivatives like esters and amides.\n'
                 '\n'
                 "Here's my improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 64,
    'num_false_positives': 100,
    'num_true_negatives': 4100,
    'num_false_negatives': 3,
    'num_negatives': None,
    'precision': 0.3902439024390244,
    'recall': 0.9552238805970149,
    'f1': 0.5541125541125541,
    'accuracy': 0.9758612608389969,
    'negative_predictive_value': 0.9992688276870583}