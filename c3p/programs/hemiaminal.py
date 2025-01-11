"""
Classifies: CHEBI:73080 hemiaminal
"""
"""
Classifies: CHEBI:60324 hemiaminal
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_hemiaminal(smiles: str):
    """
    Determines if a molecule contains a hemiaminal group based on its SMILES string.
    A hemiaminal has both an amino group and a hydroxy group attached to the same carbon atom.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a hemiaminal group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # Collection of SMARTS patterns for different hemiaminal types
    patterns = [
        # Basic hemiaminal pattern - more relaxed version
        '[CX4](-[OX2])(-[NX3])',
        
        # Cyclic hemiaminal patterns
        '[CX4;R](-[OX2])(-[NX3;R])',
        '[CX4;R](-[OX2])(-[NX3])',
        
        # Pattern for charged species
        '[CX4](-[OX2])(-[NX4+])',
        
        # Bridged bicyclic patterns common in natural products
        '[CX4;R2](-[OX2])(-[NX3;R])',
        
        # Pattern for fused ring systems
        '[CX4;R2](-[OX2;R0])(-[NX3;R])',
        
        # Pattern catching more complex cases
        '[CX4](-[OX2])(-[NX3,NX4+])',
    ]

    matches = set()
    for pattern in patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat is not None:
            these_matches = mol.GetSubstructMatches(pat)
            matches.update(match[0] for match in these_matches)  # Add carbon indices

    if not matches:
        return False, "No hemiaminal group found"

    # Patterns to exclude false positives
    exclude_patterns = [
        '[NX3](=[OX1])',  # Nitro group
        '[CX3](=O)[OX2H1]',  # Carboxylic acid
        '[NX3]-[OX2]',  # N-O bond
        '[CX3](=O)-[OX2]',  # Ester
        '[CX3]=N',  # Imine
        '[CX3](=[OX1])[NX3]'  # Amide
    ]
    
    # Remove false positives
    for pattern in exclude_patterns:
        pat = Chem.MolFromSmarts(pattern)
        if pat:
            exclude_matches = mol.GetSubstructMatches(pat)
            for match in exclude_matches:
                matches.discard(match[0])

    if not matches:
        return False, "Found potential matches but they are excluded patterns"

    # Additional validation
    valid_carbons = set()
    for carbon_idx in matches:
        carbon = mol.GetAtomWithIdx(carbon_idx)
        neighbors = carbon.GetNeighbors()
        
        # Check neighbor atoms
        n_count = 0
        o_count = 0
        for neighbor in neighbors:
            if neighbor.GetAtomicNum() == 7:  # N
                if neighbor.GetFormalCharge() >= 0:  # Allow neutral or positive N
                    n_count += 1
            elif neighbor.GetAtomicNum() == 8:  # O
                if neighbor.GetFormalCharge() == 0:  # Neutral O only
                    o_count += 1
        
        # Must have exactly one N and one O
        if n_count == 1 and o_count == 1:
            # Check hybridization of carbon
            if carbon.GetHybridization() == Chem.HybridizationType.SP3:
                valid_carbons.add(carbon_idx)

    if not valid_carbons:
        return False, "No valid hemiaminal groups found after structural validation"

    reason = f"Found {len(valid_carbons)} valid hemiaminal group(s)"
    return True, reason


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:73080',
                          'name': 'hemiaminal',
                          'definition': 'Any organic amino compound that has '
                                        'an amino group and a hydroxy group '
                                        'attached to the same carbon atom. '
                                        'Hemiaminals are intermediates in the '
                                        'formation of imines by addition of an '
                                        'amine to an aldehyde or ketone; those '
                                        'derived from primary amines are '
                                        'particularly unstable.',
                          'parents': ['CHEBI:33822', 'CHEBI:50047'],
                          'xrefs': ['Wikipedia:Hemiaminal'],
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
               '[H][C@@]12C[C@]3(CN4CC[C@@]56C(Nc7ccccc57)=C(C[C@@]5(CCO[C@@]35[H])[C@]46[H])C(=O)OC)CN3c4c(OC)cccc4[C@@]4(CCN5CC[C@]6([H])OCC[C@]6(C1)[C@@]45[H])[C@]23O '
               'NAME: vobtusine REASON: MISSED No valid hemiaminal groups '
               'found after structural validation\n'
               ' * SMILES: '
               'N1[C@@]2([C@@]3(C4=C1C=CC=C4)[C@@H](C5[C@@H]([N+]6([C@]2(C[C@H]5[C@@H]([C@H]6O)CC)[H])C)C3)O)[H] '
               'NAME: 4-methylnorajmaline REASON: MISSED No valid hemiaminal '
               'groups found after structural validation\n'
               ' * SMILES: C(CCCCCCCCCCC)(=O)N[C@H](C(O)=O)O NAME: '
               'N-dodecanoyl-(2S)-hydroxyglycine REASON: MISSED No valid '
               'hemiaminal groups found after structural validation\n'
               ' * SMILES: '
               'N12C=3C=4N(C=5C=CC=CC5C4C6=C(C3C7=C1C=CC=C7)C(NC6)=O)[C@@]8(O[C@]2(C[C@]8(N)C(=O)NC)[H])C '
               'NAME: TTT-3002 REASON: MISSED No hemiaminal group found\n'
               ' * SMILES: '
               'NC(=O)C1=CN([C@@H](O)CC1)[C@@H]1O[C@H](COP(O)(=O)OP(O)(=O)OC[C@H]2O[C@H]([C@H](O)[C@@H]2O)n2cnc3c(N)ncnc23)[C@@H](O)[C@H]1O '
               'NAME: (S)-NADHX REASON: MISSED No valid hemiaminal groups '
               'found after structural validation\n'
               ' * SMILES: '
               'C1[C@@]2(N3CC[C@@]42[C@]5(N(C=6C4=CC(=C(C6)OC)OC)C(C[C@]7([C@@]5([C@@]1(C(=CCO7)C3)[H])[H])[H])=O)[H])O '
               'NAME: pseudobrucine REASON: MISSED No valid hemiaminal groups '
               'found after structural validation\n'
               ' * SMILES: '
               '[C@@]123[C@@]4([C@]5([C@@](CC1=O)(C(CN([C@H](C2)O)C)=CCO[C@]5(CC(N4C6=C3C=CC=C6)=O)[H])[H])[H])[H] '
               'NAME: 5-hydroxyvomicine REASON: MISSED No valid hemiaminal '
               'groups found after structural validation\n'
               ' * SMILES: '
               'CN1[C@H]2CN3C4=C([C@H](COC(N)=O)[C@@]3(O)[C@@H]12)C(=O)C(O)=C(C)C4=O '
               'NAME: 7-demethylmitomycin B REASON: MISSED No valid hemiaminal '
               'groups found after structural validation\n'
               ' * SMILES: '
               '[H][C@@]12C[C@H]3[C@H](CC)[C@@H](O)N1[C@H]1C[C@@]4([C@H](OC(C)=O)C31)c1ccccc1N(C)[C@@]24[H] '
               'NAME: 17-O-acetylajmaline REASON: MISSED No valid hemiaminal '
               'groups found after structural validation\n'
               ' * SMILES: '
               'C1[C@@]2(N3CC[C@@]42[C@]5(N(C=6C4=CC=CC6)C(C[C@]7([C@@]5([C@@]1(C(=CCO7)C3)[H])[H])[H])=O)[H])O '
               'NAME: pseudostrychnine REASON: MISSED No valid hemiaminal '
               'groups found after structural validation\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'OC(=O)CCNC(O)=O',
                                     'name': 'N-carboxy-beta-alanine',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': 'C[C@@H]1CN([C@@H](COC2=C(C=CC(=C2)NC(=O)NC3=CC=CC=C3F)C(=O)N(C[C@H]1OC)C)C)CCC(F)(F)F',
                                     'name': '1-(2-fluorophenyl)-3-[(5S,6R,9R)-5-methoxy-3,6,9-trimethyl-2-oxo-8-(3,3,3-trifluoropropyl)-11-oxa-3,8-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]urea',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': 'O(C1=CC(=C(O)C=C1)C2=C(O)C=CC(=C2)OC)C',
                                     'name': '5,5′-dimethoxybiphenyl-2,2′-diol',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': 'COC1=CC=C(CCNC[C@H](O)C2=CC=C(O)C=C2)C=C1OC',
                                     'name': 'denopamine',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': '[Zn++].[S-]C(=S)NCCNC([S-])=S',
                                     'name': 'zineb',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': 'P(OC[C@H]1O[C@@H](N2C=CC(=NC2=O)N)C(O)[C@H]1O)(OP(OC[C@H](OC(=O)CCCCCCCCCCCC(C)C)COC(=O)CCCCCCCCC(CC)C)(O)=O)(O)=O',
                                     'name': 'CDP-DG(a-13:0/i-15:0)',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': 'O=C1OC(O)C2=C1C[C@](O)([C@H]3CC(C[C@H]3[C@@H]2O)(C)C)C',
                                     'name': 'Lactarolide A',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': 'O([C@H]1[C@H](O)[C@H](OC(O)[C@@H]1NC(=O)C)CO[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO)[C@@H]3O[C@@H]([C@H](O)[C@H](O)[C@H]3O)CO',
                                     'name': 'N-[(3R,4R,5S,6R)-2,5-Dihydroxy-4-[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-6-[[(2R,3R,4S,5R,6R)-3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-3-yl]acetamide',
                                     'reason': 'No hemiaminal group found'},
                                 {   'smiles': 'COc1ccc(\\C=C/C2CCC=CC2c2ccc(OC)c(OC)c2)cc1OC',
                                     'name': '3-(3,4-dimethoxyphenyl)-4-[(Z)-3,4-dimethoxystyryl]cyclohex-1-ene',
                                     'reason': 'No hemiaminal group found'}],
    'sample_false_negatives': [   {   'smiles': 'N12C=3C=4N(C=5C=CC=CC5C4C6=C(C3C7=C1C=CC=C7)C(NC6)=O)[C@@]8(O[C@]2(C[C@]8(N)C(=O)NC)[H])C',
                                      'name': 'TTT-3002',
                                      'reason': 'No hemiaminal group found'},
                                  {   'smiles': 'CCCOc1ccc2n3[C@H]4C[C@](O)(C(=O)OC)[C@](C)(O4)n4c5ccccc5c5c6CNC(=O)c6c(c2c1)c3c45',
                                      'name': 'KT 5926',
                                      'reason': 'No hemiaminal group found'},
                                  {   'smiles': 'COC(=O)[C@]1(C[C@H]2O[C@]1(C)n1c3ccccc3c3c4CN(C)C(=O)c4c4c5ccccc5n2c4c13)OC',
                                      'name': 'KT 5823',
                                      'reason': 'No hemiaminal group found'},
                                  {   'smiles': '[H][C@]12N3CCC[C@@]1(CC)C[C@](O)(C(=O)OC)n1c2c(CC3)c2ccccc12',
                                      'name': 'vincamine',
                                      'reason': 'No hemiaminal group found'}],
    'attempt': 3,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. Pattern Matching Issues:\n'
                 '- The previous code was too restrictive in its SMARTS '
                 'patterns\n'
                 '- It required explicit hydrogens on the hydroxyl group '
                 '(-[OX2H1]) which might miss some cases\n'
                 '- The nitrogen pattern was too general and might need more '
                 'specific valence states\n'
                 '\n'
                 '2. False Negatives Analysis:\n'
                 'Looking at the failed examples:\n'
                 '- Many are complex natural products with bridged structures\n'
                 '- Some have the hemiaminal as part of a ring system\n'
                 '- The hydroxy group may not always have explicit hydrogens\n'
                 '- Some structures have the hemiaminal nitrogen as part of a '
                 'larger ring system\n'
                 '\n'
                 '3. Improvements Needed:\n'
                 '- Relax the hydroxyl pattern to catch more cases\n'
                 '- Add patterns for ring-based hemiaminals\n'
                 '- Better handle bridged bicyclic systems\n'
                 '- Consider charged species (like in 4-methylnorajmaline)\n'
                 '- Add specific patterns for natural product scaffolds\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 22,
    'num_false_positives': 100,
    'num_true_negatives': 8849,
    'num_false_negatives': 4,
    'num_negatives': None,
    'precision': 0.18032786885245902,
    'recall': 0.8461538461538461,
    'f1': 0.29729729729729726,
    'accuracy': 0.9884122562674095,
    'negative_predictive_value': 0.9995481757596295}