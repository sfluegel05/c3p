"""
Classifies: CHEBI:74716 withanolide
"""
"""
Classifies: withanolide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_withanolide(smiles: str):
    """
    Determines if a molecule is a withanolide based on its SMILES string.
    A withanolide is defined as any steroid lactone that is a C28 steroid 
    with a modified side chain forming a lactone ring and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a withanolide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for tetracyclic steroid nucleus (fused ring system with 4 rings: 6-6-6-5)
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    bond_rings = ring_info.BondRings()
    
    # Build a list of rings and their sizes
    ring_sizes = [len(ring) for ring in atom_rings]
    num_rings = len(atom_rings)
    
    # Build a connectivity map of rings
    ring_adj_matrix = [[False]*num_rings for _ in range(num_rings)]
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            # If rings share bonds, they are fused
            if set(bond_rings[i]) & set(bond_rings[j]):
                ring_adj_matrix[i][j] = True
                ring_adj_matrix[j][i] = True

    # Function to find connected components (fused ring systems)
    def get_fused_systems():
        visited = set()
        fused_systems = []
        for i in range(num_rings):
            if i not in visited:
                stack = [i]
                fused_system = []
                while stack:
                    idx = stack.pop()
                    if idx not in visited:
                        visited.add(idx)
                        fused_system.append(idx)
                        for j in range(num_rings):
                            if ring_adj_matrix[idx][j] and j not in visited:
                                stack.append(j)
                fused_systems.append(fused_system)
        return fused_systems

    fused_systems = get_fused_systems()
    has_steroid_nucleus = False
    steroid_atoms = set()

    for system in fused_systems:
        if len(system) == 4:
            system_ring_sizes = [ring_sizes[idx] for idx in system]
            system_ring_sizes.sort()
            if system_ring_sizes == [5,6,6,6]:
                has_steroid_nucleus = True
                # Collect atoms in the steroid nucleus
                for idx in system:
                    steroid_atoms.update(atom_rings[idx])
                break

    if not has_steroid_nucleus:
        return False, "No tetracyclic steroid nucleus (rings of sizes 6-6-6-5) found"

    # Check for lactone ring in side chain (cyclic ester not part of steroid nucleus)
    # Define lactone pattern (cyclic ester)
    lactone_pattern = Chem.MolFromSmarts('C(=O)O[C;R]')
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    has_lactone = False

    for match in lactone_matches:
        lactone_o = match[2]
        # Check if the oxygen is part of a ring
        for ring in atom_rings:
            if lactone_o in ring:
                # Check if ring is not part of steroid nucleus
                ring_atoms_set = set(ring)
                if not ring_atoms_set.issubset(steroid_atoms):
                    has_lactone = True
                    break
        if has_lactone:
            break

    if not has_lactone:
        return False, "No lactone ring found in side chain"

    # Check number of carbons (should be around 28)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 26 or c_count > 30:
        return False, f"Carbon count ({c_count}) not consistent with C28 steroid"

    return True, "Molecule is a withanolide (tetracyclic steroid nucleus with side-chain lactone ring)"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:74716',
                          'name': 'withanolide',
                          'definition': 'Any steroid lactone that is a C28 '
                                        'steroid with a modified side chain '
                                        'forming a lactone ring and its '
                                        'substituted derivatives.',
                          'parents': ['CHEBI:26766'],
                          'xrefs': [   'LIPID_MAPS_instance:LMST01160000',
                                       'Wikipedia:Withanolide'],
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
               'CC1=C(C(=O)O[C@@H](C1)C(C)[C@H]2CC[C@H]3[C@]2(CCC4[C@@H]3C[C@@H]5[C@@]6([C@]4(C(=O)C=C[C@@H]6O)C)O5)C)CO '
               'NAME: LSM-5073 REASON: MISSED No steroid nucleus found\n'
               ' * SMILES: '
               'O=C1C2=C(C=C[C@]3([C@]2(CC[C@@H]3[C@@H](C[C@H]4OC(=O)C(=C4C)C)C)C)C)[C@]5(CC[C@H](C([C@@H]5C1)(C)C)OC(=O)C[C@@](O)(CC(=O)OC)C)C '
               'NAME: Palustrisolide G REASON: MISSED No steroid nucleus '
               'found\n'
               ' * SMILES: '
               'C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)[C@H](CC(=O)[C@]4(C)[C@H]3CC[C@]12C)OS(O)(=O)=O)[C@H]1CC(C)=C(CO)C(=O)O1 '
               'NAME: 2,3-dihydro-3beta-O-sulfate withaferin A REASON: MISSED '
               'No steroid nucleus found\n'
               ' * SMILES: '
               'O=C1O[C@]2(OC[C@@]34[C@](C5=C([C@@]6([C@H](C([C@H](O)CC6)(C)C)CC5)C)C[C@@H]3O)(CC[C@@H]4[C@@H](C2)C)C)[C@H]([C@@H]1C)C '
               'NAME: Fomefficinol A REASON: MISSED No steroid nucleus found\n'
               ' * SMILES: '
               'ClC1CC2C3(O)C(C(O)(CC3)C(O)(C4OC(=O)C(=C(C4)C)C)C)(CCC2C5(C1(O)C(O)C=CC5=O)C)C '
               'NAME: physalolactone REASON: MISSED No steroid nucleus found\n'
               ' * SMILES: '
               'O=C1C=C2O[C@@H]3C[C@@]4([C@H]([C@@H](OC(=O)C)C[C@]4(C5=C3[C@]2(C)[C@H](C1(C)C)CC5)C)[C@@H](C[C@@H]6OC(=O)C(=C6)C)C)C '
               'NAME: Lepiotaprocerin E REASON: MISSED No steroid nucleus '
               'found\n'
               ' * SMILES: '
               'OC12C(C(CC1)C(O)(C3OC(=O)C(=C(C3)C)CO)C)(CCC4C2CC(O)C5(O)C4(C(=O)CCC5)C)C '
               'NAME: '
               '(5alpha,6beta,14alpha,20R,22R)-5,6,14,20,27-Pentahydroxy-1-oxowith-24-enolide '
               'REASON: MISSED No steroid nucleus found\n'
               ' * SMILES: '
               'O=C1C2=C([C@@]3([C@@H](OC(=O)C)C[C@]4([C@]3(C1)C)O[C@]5(OC(=O)[C@H](C5)C)C[C@H]4C)C)CC[C@@H]6[C@@]2(CCC(C6(C)C)=O)C '
               'NAME: Spiroganocalitone C REASON: MISSED No steroid nucleus '
               'found\n'
               ' * SMILES: '
               'O=C1OC([C@H]([C@@H]2[C@@]3([C@@](C4=C([C@@]5([C@H](C([C@@H](O)CC5)(C)C)CC4)COC(=O)C)CC3)(C)CC2)C)C)CC=C1C '
               'NAME: Colobetaolactone B REASON: MISSED No steroid nucleus '
               'found\n'
               ' * SMILES: '
               'O1C2(C3(C(O)(C(OC(=O)C)C12)C4C(CC3)C5(C(O)(C(O)C4)CC=CC5=O)C)C)C(C6OC(=O)C(=C(C6)C)C)C '
               'NAME: Physagulin F REASON: MISSED No steroid nucleus found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'CNC(=O)CC[C@H](N)C(O)=O',
                                     'name': 'N(5)-methyl-L-glutamine',
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'},
                                 {   'smiles': 'C[C@H](OP(O)(=O)OP(O)(=O)OC[C@H]1O[C@H]([C@H](O)[C@@H]1O)n1cnc2c1nc(N)[nH]c2=O)C(O)=O',
                                     'name': "L-lactyl-2-diphospho-5'-guanosine",
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'},
                                 {   'smiles': '[H][C@@]1(CC[C@]2(C)[C@]1([H])CC[C@]1([H])[C@@]3(C)CCCC(C)(C)[C@]3([H])CC[C@@]21C)C(=C)CCC=C(C)C',
                                     'name': 'dammara-20,24-diene',
                                     'reason': 'No lactone ring found in side '
                                               'chain'},
                                 {   'smiles': 'O=C1C2=C(O)C3=C(O)C=CC=C3C(=C2C(C(=O)C)C(C1)(O)C)C=4C=5C(C(O)=C6C4C(C(=O)C)C(O)(C)CC6=O)=C(O)C=CC5',
                                     'name': 'A 39183A',
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'},
                                 {   'smiles': 'C[C@H](O)[C@@H](O)C1=Nc2c(NC1)[nH]c(N)nc2=O',
                                     'name': 'L-threo-7,8-dihydrobiopterin',
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CN([C@H](COC2=C(C=C(C=C2)NC(=O)C)C(=O)N(C[C@H]1OC)C)C)CC3=CC=C(C=C3)F',
                                     'name': 'N-[(4S,7R,8S)-5-[(4-fluorophenyl)methyl]-8-methoxy-4,7,10-trimethyl-11-oxo-2-oxa-5,10-diazabicyclo[10.4.0]hexadeca-1(12),13,15-trien-14-yl]acetamide',
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'},
                                 {   'smiles': 'C(CCN1CCCCC1)(CCN(C(C)C)C(C)=O)(C(N)=O)C2=C(Cl)C=CC=C2',
                                     'name': 'bidisomide',
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'},
                                 {   'smiles': 'OCCCCCCCCCCC#CC=C',
                                     'name': '13-Tetradece-11-yn-1-ol',
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'},
                                 {   'smiles': 'C[C@@H]1CCCCO[C@@H]([C@@H](CN(C(=O)C2=C(O1)C=CC(=C2)N(C)C)[C@@H](C)CO)C)CN(C)CC3=CC=C(C=C3)OC',
                                     'name': '(3R,9S,10R)-16-(dimethylamino)-12-[(2S)-1-hydroxypropan-2-yl]-9-[[(4-methoxyphenyl)methyl-methylamino]methyl]-3,10-dimethyl-2,8-dioxa-12-azabicyclo[12.4.0]octadeca-1(14),15,17-trien-13-one',
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'},
                                 {   'smiles': 'O(CCC12CC3CC(C1)CC(C2)C3)C(=O)CC4=CC=C(OCC(O)CNC(C)C)C=C4',
                                     'name': 'Adaprolol',
                                     'reason': 'No tetracyclic steroid nucleus '
                                               '(rings of sizes 6-6-6-5) '
                                               'found'}],
    'sample_false_negatives': [   {   'smiles': 'CC1=C(C(=O)O[C@@H](C1)C(C)[C@H]2CC[C@H]3[C@]2(CCC4[C@@H]3C[C@@H]5[C@@]6([C@]4(C(=O)C=C[C@@H]6O)C)O5)C)CO',
                                      'name': 'LSM-5073',
                                      'reason': 'No tetracyclic steroid '
                                                'nucleus (rings of sizes '
                                                '6-6-6-5) found'},
                                  {   'smiles': 'O=C1C2=C(C=C[C@]3([C@]2(CC[C@@H]3[C@@H](C[C@H]4OC(=O)C(=C4C)C)C)C)C)[C@]5(CC[C@H](C([C@@H]5C1)(C)C)OC(=O)C[C@@](O)(CC(=O)OC)C)C',
                                      'name': 'Palustrisolide G',
                                      'reason': 'Carbon count (38) not '
                                                'consistent with C28 steroid'},
                                  {   'smiles': 'C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)[C@H](CC(=O)[C@]4(C)[C@H]3CC[C@]12C)OS(O)(=O)=O)[C@H]1CC(C)=C(CO)C(=O)O1',
                                      'name': '2,3-dihydro-3beta-O-sulfate '
                                              'withaferin A',
                                      'reason': 'No tetracyclic steroid '
                                                'nucleus (rings of sizes '
                                                '6-6-6-5) found'},
                                  {   'smiles': 'O=C1O[C@]2(OC[C@@]34[C@](C5=C([C@@]6([C@H](C([C@H](O)CC6)(C)C)CC5)C)C[C@@H]3O)(CC[C@@H]4[C@@H](C2)C)C)[C@H]([C@@H]1C)C',
                                      'name': 'Fomefficinol A',
                                      'reason': 'No tetracyclic steroid '
                                                'nucleus (rings of sizes '
                                                '6-6-6-5) found'},
                                  {   'smiles': 'O=C1C=C2O[C@@H]3C[C@@]4([C@H]([C@@H](OC(=O)C)C[C@]4(C5=C3[C@]2(C)[C@H](C1(C)C)CC5)C)[C@@H](C[C@@H]6OC(=O)C(=C6)C)C)C',
                                      'name': 'Lepiotaprocerin E',
                                      'reason': 'No tetracyclic steroid '
                                                'nucleus (rings of sizes '
                                                '6-6-6-5) found'},
                                  {   'smiles': 'O=C1C2=C([C@@]3([C@@H](OC(=O)C)C[C@]4([C@]3(C1)C)O[C@]5(OC(=O)[C@H](C5)C)C[C@H]4C)C)CC[C@@H]6[C@@]2(CCC(C6(C)C)=O)C',
                                      'name': 'Spiroganocalitone C',
                                      'reason': 'Carbon count (32) not '
                                                'consistent with C28 steroid'},
                                  {   'smiles': 'O=C1OC([C@H]([C@@H]2[C@@]3([C@@](C4=C([C@@]5([C@H](C([C@@H](O)CC5)(C)C)CC4)COC(=O)C)CC3)(C)CC2)C)C)CC=C1C',
                                      'name': 'Colobetaolactone B',
                                      'reason': 'Carbon count (32) not '
                                                'consistent with C28 steroid'},
                                  {   'smiles': 'O1C2(C3(C(O)(C(OC(=O)C)C12)C4C(CC3)C5(C(O)(C(O)C4)CC=CC5=O)C)C)C(C6OC(=O)C(=C(C6)C)C)C',
                                      'name': 'Physagulin F',
                                      'reason': 'No tetracyclic steroid '
                                                'nucleus (rings of sizes '
                                                '6-6-6-5) found'},
                                  {   'smiles': 'OC(C1C2(C(C3C(C4(C(=CC3)CC(OC5OC(C(O)C(O)C5O)CO)CC4OC(=O)C)C)CC2)CC1)C)(C6OC(=O)C(=C(C6)C)C)C',
                                      'name': 'Physalolactone B 3-glucoside',
                                      'reason': 'Carbon count (36) not '
                                                'consistent with C28 steroid'},
                                  {   'smiles': 'O1C23C(C4C(C5(O)C(C(O)(CC5)C(O)(C6OC(=O)C(=C(C6)C)C)C)(CC4)C)CC12)(C)C(=O)C=CC3=O',
                                      'name': 'Withaperuvin E',
                                      'reason': 'No tetracyclic steroid '
                                                'nucleus (rings of sizes '
                                                '6-6-6-5) found'}],
    'attempt': 3,
    'reasoning': "In analyzing the previous code and the given outcomes, it's "
                 'clear that the function `is_withanolide` failed to recognize '
                 'the steroid nucleus in any of the provided SMILES strings, '
                 'consistently returning "No steroid nucleus found".\n'
                 '\n'
                 '**Possible Reasons for Failure:**\n'
                 '\n'
                 '1. **Overly Specific Steroid Pattern:**\n'
                 '   - The SMARTS pattern used for the steroid nucleus '
                 "(`'C1CCC2C(C1)CCC3C2CCC4C3(CCCC4)C'`) is too specific.\n"
                 '   - This pattern represents a specific steroid skeleton and '
                 'does not account for variations commonly found in natural '
                 'withanolides, such as double bonds, ring fusions, '
                 'stereochemistry, and various substituents.\n'
                 '\n'
                 '2. **Lactone Ring Pattern Limitations:**\n'
                 "   - The lactone pattern (`'C1OC(=O)[C;R]1'`) focuses only "
                 'on 5-membered lactone rings.\n'
                 '   - Withanolides may have lactone rings of different sizes '
                 'and complexities in their side chains.\n'
                 '\n'
                 '**Improvements:**\n'
                 '\n'
                 '1. **Generalize Steroid Nucleus Detection:**\n'
                 '   - Use a more flexible approach to identify the '
                 'tetracyclic steroid nucleus.\n'
                 '   - Instead of a rigid SMARTS pattern, analyze the '
                 "molecule's ring structures to identify a fused ring system "
                 'with four rings, typically of sizes 6-6-6-5.\n'
                 '\n'
                 '2. **Generalize Lactone Ring Detection:**\n'
                 '   - Use a more general method to detect lactone rings, '
                 'regardless of ring size.\n'
                 '   - Ensure the lactone ring is part of the side chain and '
                 'not within the steroid nucleus.\n'
                 '\n'
                 '3. **Analyze Ring Connectivity:**\n'
                 "   - Utilize RDKit's ring information to identify fused ring "
                 'systems and their sizes.\n'
                 '   - This allows for more flexibility in detecting the core '
                 'structure of withanolides.\n'
                 '\n'
                 '4. **Consider Carbon Count:**\n'
                 '   - Verify that the molecule has approximately 28 carbons, '
                 'as withanolides are C28 steroids.\n'
                 '\n'
                 '**Revised Code:**\n'
                 '\n'
                 'Below is the improved code that implements these changes:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 32,
    'num_false_positives': 71,
    'num_true_negatives': 142103,
    'num_false_negatives': 94,
    'num_negatives': None,
    'precision': 0.3106796116504854,
    'recall': 0.25396825396825395,
    'f1': 0.2794759825327511,
    'accuracy': 0.9988404778636683,
    'negative_predictive_value': 0.9993389452660746}