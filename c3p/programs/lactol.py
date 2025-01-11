"""
Classifies: CHEBI:38131 lactol
"""
"""
Classifies: CHEBI:35879 lactol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_lactol(smiles: str):
    """
    Determines if a molecule is a lactol based on its SMILES string.
    Lactols are cyclic hemiacetals formed by intramolecular addition of a hydroxy 
    group to an aldehydic or ketonic carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a lactol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic cyclic hemiacetal patterns
    patterns = [
        # Ring oxygen connected to carbon with OH or OR
        "[O;R][C;R]([OH1,OR])",
        
        # Ring oxygen connected to carbon with OH, more specific
        "[O;R][C;R]([OH1])",
        
        # Hemiacetal in equilibrium with ketone/aldehyde form
        "[O;R][C;R](=O)",
        
        # Common natural product lactol pattern
        "[O;R][C;R]([OH1,OR])[C,H]",
        
        # Bridged lactol pattern
        "[O;R][C;R]([OH1,OR])[C;R]",
        
        # Fused ring lactol pattern
        "[O;R][C;R]([OH1,OR])[C;R0]"
    ]

    # Convert patterns to RDKit SMARTS objects
    smarts_patterns = [Chem.MolFromSmarts(p) for p in patterns]

    # Check for matches
    matches = []
    for pattern in smarts_patterns:
        if pattern is not None:  # Ensure valid SMARTS pattern
            matches.extend(mol.GetSubstructMatches(pattern))

    if not matches:
        return False, "No lactol structure found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    
    # Verify that matched structures are part of rings
    valid_lactol = False
    for match in matches:
        # Check first two atoms (O and C) are in same ring
        rings_with_o = set(i for i, ring in enumerate(ring_info.AtomRings()) 
                          if match[0] in ring)
        rings_with_c = set(i for i, ring in enumerate(ring_info.AtomRings()) 
                          if match[1] in ring)
        
        if rings_with_o.intersection(rings_with_c):
            valid_lactol = True
            break

    if not valid_lactol:
        return False, "Matched atoms not in same ring"

    # Additional checks for specific cases
    # Pattern to identify simple glycosides
    glycoside = Chem.MolFromSmarts("[OR0][C;R]1[O;R][C;R][C;R][C;R][C;R]1")
    
    # If molecule only matches glycoside pattern and has typical sugar structure,
    # it's probably a glycoside rather than a lactol
    if (mol.HasSubstructMatch(glycoside) and 
        len(matches) == 1 and 
        any(a.GetSymbol() == 'O' for a in mol.GetAtoms()) > 5):
        return False, "Structure appears to be a glycoside rather than a lactol"

    return True, "Contains cyclic hemiacetal structure characteristic of lactols"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:38131',
                          'name': 'lactol',
                          'definition': 'Cyclic hemiacetals formed by '
                                        'intramolecular addition of a hydroxy '
                                        'group to an aldehydic or ketonic '
                                        'carbonyl group. They are thus '
                                        '1-oxacycloalkan-2-ols or unsaturated '
                                        'analogues.',
                          'parents': ['CHEBI:5653'],
                          'xrefs': ['Wikipedia:Lactol'],
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
               'O1C2=C(C(CC1(C3=CC=CC=C3)O)=O)C(=CC(=C2)O)O NAME: '
               '2,5,7-trihydroxy-2-phenyl-2,3-dihydro-4H-chromen-4-one REASON: '
               'MISSED No lactol structure found\n'
               ' * SMILES: '
               'OCC[C@@H]1[C@@H](O)Oc2cc3C(=O)c4cc(O)cc(O)c4C(=O)c3c(O)c12 '
               'NAME: (2S-3S)-versiconal hemiacetal REASON: MISSED No lactol '
               'structure found\n'
               ' * SMILES: '
               'CO[C@H]1\\C=C\\C=C(C)\\C[C@H](C)[C@H](O)[C@H](C)\\C=C(/C)\\C=C(C)\\C(=O)O[C@@H]1[C@@H](C)[C@@H](O)[C@H](C)[C@@]1(O)C[C@@H](OC(=O)\\C=C\\C(O)=O)[C@@H](C)[C@@H](C)O1 '
               'NAME: hygrolidin REASON: MISSED No lactol structure found\n'
               ' * SMILES: '
               'C/1(\\C[C@H]2O[C@@H](C1)C[C@]3(C([C@H](C[C@@H](C[C@H](CC(O[C@@H]([C@@H](C)O)C[C@@]4(C\\C(\\[C@@H]([C@](C(C=C2)(C)C)(O4)O)OC(/C=C/C=C/CCC)=O)=C/C(OC)=O)[H])=O)O)O3)O)(C)C)O)=C\\C(OC)=O '
               'NAME: bryostatin 2 REASON: MISSED No lactol structure found\n'
               ' * SMILES: '
               'CC[C@@H]1O[C@H](CC[C@@]1(C)O)[C@]1(C)C[C@@H](C)[C@@](O)(O1)[C@]1(C)CC[C@H](O1)[C@@H](O)[C@]1(O)O[C@H](C[C@@]2(O)O[C@H](CC[C@@H]2C)[C@@H](C)[C@H](O)C[C@@H]2O[C@H]([C@@H](C)C(O)=O)[C@@H](C)C[C@H]2C)[C@@H](C)C[C@H]1C '
               'NAME: alborixin REASON: MISSED No lactol structure found\n'
               ' * SMILES: '
               'CCC(O)CC1CCCC2(CC3OC(=O)\\C=C\\C(C)(O)C(O)C(C)C(O)C(OC4CCC(C(C)O4)N(C)C)C(O)C(C)(O)CCCCC\\C=C\\C4CC(C)(C)OC4(O)CC(O2)C3C)O1 '
               'NAME: NK154183B REASON: MISSED No lactol structure found\n'
               ' * SMILES: '
               'CC(=O)OCC[C@@H]1[C@@H](O)Oc2cc3C(=O)c4cc(O)cc(O)c4C(=O)c3c(O)c12 '
               'NAME: (2S,3S)-versiconal hemiacetal acetate REASON: MISSED No '
               'lactol structure found\n'
               ' * SMILES: '
               'C/1(\\C[C@H]2O[C@@H](C1)C[C@]3(C([C@H](C[C@@H](C[C@H](CC(O[C@@H]([C@@H](C)O)C[C@@]4(C\\C(\\[C@@H]([C@](C(C=C2)(C)C)(O4)O)OC(/C=C/C=C/CCC)=O)=C/C(OC)=O)[H])=O)O)O3)OC(C)=O)(C)C)O)=C\\C(OC)=O '
               'NAME: bryostatin 1 REASON: MISSED No lactol structure found\n'
               ' * SMILES: '
               '[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@@]23[C@]4([C@@](C2)([C@@]5(O[C@]4(O[C@]3(C5)C)[H])O)[H])COC(C6=CC=C(C=C6)O)=O '
               'NAME: oxypaeoniflorin REASON: MISSED No lactol structure '
               'found\n'
               ' * SMILES: OC[C@]1(O)OC[C@@H](O)[C@@H](O)[C@@H]1O NAME: '
               'alpha-D-fructopyranose REASON: MISSED No lactol structure '
               'found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CC(=O)O[C@H]1CC[C@]23C[C@@H]1OO[C@@]2(C)C(=O)CCC3(C)C',
                                     'name': 'Talaperoxide B',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'C1=CC=C(C(=C1)C#CC2=CC=C(C=C2)[C@H]3[C@H](N[C@@H]3C#N)CO)F',
                                     'name': '(2S,3R,4S)-3-[4-[2-(2-fluorophenyl)ethynyl]phenyl]-4-(hydroxymethyl)-2-azetidinecarbonitrile',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'OC(=O)Cc1cn(nc1-c1ccc(Cl)cc1)-c1ccccc1',
                                     'name': 'lonazolac',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'C(CCCCCCCC(=O)O)CCCCCCCCC(=O)O',
                                     'name': 'octadecanedioic acid',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'C\\C=C(/C)CC\\C=C(/C)C(=O)OCC(C)(C)CC1=C(O)C(=O)c2ccccc2C1=O',
                                     'name': 'rhinacanthin C',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'COC1=C(C=C2C(=C1)C(=NC=N2)NC3=CC(=C(C(=C3)Br)O)Br)OC',
                                     'name': '2,6-dibromo-4-[(6,7-dimethoxy-4-quinazolinyl)amino]phenol',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'S(=O)(CC1=CC=CC=C1)C',
                                     'name': 'Methyl benzyl sulfoxide',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'C=CCOC1=NS(=O)(=O)c2ccccc12',
                                     'name': 'probenazole',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'C(=O)([C@@H](NC(=O)CC)CCSC)[O-]',
                                     'name': 'N-propanoyl-L-methioninate',
                                     'reason': 'No lactol structure found'},
                                 {   'smiles': 'O1C(C(O)C(O)C(O)C1OC(CCC=2C(CC(O)CC2C)(C)C)C)COC3OCC(O)(C3O)CO',
                                     'name': '(3b,9R)-5-Megastigmene-3,9-diol '
                                             '9-[apiosyl-(1->6)-glucoside]',
                                     'reason': 'No lactol structure found'}],
    'sample_false_negatives': [   {   'smiles': 'C=12C(=C3C(=CC(=CC3=CC1C=C(OC2=O)CC(OC)=O)OC)O)O',
                                      'name': 'methyl '
                                              '2-[(3S)-9,10-dihydroxy-7-methoxy-1-oxo-1H,3H,4H-naphtho[2,3-c]pyran-3-yl]acetate',
                                      'reason': 'No lactol structure found'},
                                  {   'smiles': 'C=12C(=C3C(=CC(=CC3=CC1C=C(OC2=O)CC(C)=O)OC)O)O',
                                      'name': '9,10-dihydroxy-7-methoxy-3-(2-oxopropyl)-1H-benzo[g]isochromen-1-one',
                                      'reason': 'No lactol structure found'}],
    'attempt': 2,
    'reasoning': 'Let me analyze the previous attempt and suggest '
                 'improvements:\n'
                 '\n'
                 '1. The main issue seems to be that the SMARTS patterns were '
                 'too restrictive and missed many valid lactol structures. The '
                 'patterns were looking for very specific arrangements:\n'
                 '   - `[O;R{5-6}][C;R{5-6}][OH1]` - requiring a hydroxyl '
                 'directly on the ring carbon\n'
                 '   - `[O;R{5-6}][C;R{5-6}][O;R{5-6}]` - requiring both '
                 'oxygens to be in the ring\n'
                 '   - `[O;R{5-6}][C;R{5-6}](=O)` - looking for ring carbonyl\n'
                 '\n'
                 '2. Looking at the failed examples:\n'
                 '   - Many contain more complex ring systems where the lactol '
                 'is part of a larger structure\n'
                 '   - Some have substituents or variations in the ring that '
                 "weren't captured\n"
                 '   - The ring size restriction to 5-6 members may have been '
                 'too limiting\n'
                 '\n'
                 '3. Improvements needed:\n'
                 '   - Broaden the SMARTS patterns to capture more variations\n'
                 '   - Include patterns for bridged and fused ring systems\n'
                 '   - Consider the hemiacetal nature more generally (O-C-OH '
                 'or O-C-OR where C is part of ring)\n'
                 '   - Remove strict ring size limitations\n'
                 '   - Add patterns for common lactol-containing natural '
                 'products\n'
                 '\n'
                 "Here's my improved version:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 101,
    'num_false_positives': 100,
    'num_true_negatives': 496,
    'num_false_negatives': 2,
    'num_negatives': None,
    'precision': 0.5024875621890548,
    'recall': 0.9805825242718447,
    'f1': 0.6644736842105263,
    'accuracy': 0.8540772532188842,
    'negative_predictive_value': 0.9959839357429718}