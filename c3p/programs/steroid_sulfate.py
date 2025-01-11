"""
Classifies: CHEBI:16158 steroid sulfate
"""
"""
Classifies: steroid sulfate
"""

from rdkit import Chem

def get_fused_ring_subgraphs(mol):
    """
    Helper function to identify fused ring systems in a molecule.

    Args:
        mol: rdkit.Chem.Mol object

    Returns:
        List of sets, each set contains indices of rings that are fused together
    """
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Build adjacency lists for rings
    ring_adj_list = []
    num_rings = len(atom_rings)
    for i in range(num_rings):
        ring_i_atoms = set(atom_rings[i])
        fused = []
        for j in range(num_rings):
            if i != j:
                ring_j_atoms = set(atom_rings[j])
                if ring_i_atoms & ring_j_atoms:
                    fused.append(j)
        ring_adj_list.append(fused)

    # Find connected components in ring adjacency graph
    from collections import deque
    visited = set()
    fused_ring_systems = []
    for i in range(num_rings):
        if i not in visited:
            system = set()
            queue = deque([i])
            while queue:
                idx = queue.popleft()
                if idx not in visited:
                    visited.add(idx)
                    system.add(idx)
                    neighbors = ring_adj_list[idx]
                    for n in neighbors:
                        if n not in visited:
                            queue.append(n)
            fused_ring_systems.append(system)
    return fused_ring_systems

def is_steroid_sulfate(smiles: str):
    """
    Determines if a molecule is a steroid sulfate based on its SMILES string.
    A steroid sulfate is a steroid molecule where at least one hydroxy group is esterified with sulfuric acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid sulfate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Get fused ring systems
    fused_ring_systems = get_fused_ring_subgraphs(mol)

    # Identify steroid backbone
    steroid_found = False
    steroid_atom_indices = set()

    for ring_idxs in fused_ring_systems:
        # Collect atom indices and ring sizes
        ring_atoms = set()
        ring_sizes = []
        for ring_idx in ring_idxs:
            ring = atom_rings[ring_idx]
            ring_atoms.update(ring)
            ring_sizes.append(len(ring))
        # Check for 4-ring system with sizes 6,6,6,5 and 17 carbons
        if len(ring_idxs) == 4:
            ring_sizes_sorted = sorted(ring_sizes)
            if ring_sizes_sorted == [5,6,6,6] and len(ring_atoms) == 17:
                steroid_found = True
                steroid_atom_indices = ring_atoms
                break

    if not steroid_found:
        return False, "No steroid backbone found"

    # Define sulfate ester group patterns (protonated and deprotonated)
    sulfate_pattern = Chem.MolFromSmarts('OS(=O)(=O)[O-]')
    sulfate_pattern_unprotonated = Chem.MolFromSmarts('OS(=O)(=O)O')

    sulfate_matches = mol.GetSubstructMatches(sulfate_pattern) + mol.GetSubstructMatches(sulfate_pattern_unprotonated)

    if not sulfate_matches:
        return False, "No sulfate ester group found"

    # Check if sulfate is attached to steroid backbone via oxygen
    for match in sulfate_matches:
        ester_oxygen_idx = match[0]  # Oxygen connected to sulfur
        ester_oxygen_atom = mol.GetAtomWithIdx(ester_oxygen_idx)
        # Check neighbors of ester oxygen atom
        for neighbor in ester_oxygen_atom.GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx in steroid_atom_indices:
                return True, "Contains steroid backbone with sulfate ester group attached via oxygen"

    return False, "Sulfate group not attached to steroid backbone via ester linkage"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16158',
                          'name': 'steroid sulfate',
                          'definition': 'A sulfuric ester obtained by the '
                                        'formal condensation of a hydroxy '
                                        'group of any steroid with sulfuric '
                                        'acid.',
                          'parents': [   'CHEBI:25704',
                                         'CHEBI:26819',
                                         'CHEBI:47880'],
                          'xrefs': ['KEGG:C02590'],
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
               'S(O[C@@H]1C[C@]2([C@]3([C@@]([C@@]4(C(CC3)CC(=O)C=C4O)C)(CC[C@@]2(C1)C)[H])[H])[H])(O)(=O)=O '
               'NAME: 3b,16a-Dihydroxyandrostenone sulfate REASON: MISSED No '
               'steroid backbone found\n'
               ' * SMILES: '
               'C1[C@@]2([C@@]([C@@]3(C(C[C@@H](OS(O)(=O)=O)CC3)=C1)C)(CC[C@]4([C@]2(C[C@H](C4=O)O)[H])C)[H])[H] '
               'NAME: 16alpha-hydroxydehydroepiandrosterone 3-sulfate REASON: '
               'MISSED No steroid backbone found\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@]2([H])[C@]1(C)CC[C@]1([H])[C@@]34CC[C@@H](C[C@]3(O)[C@H](O)C[C@@]21OC4)OS(O)(=O)=O)[C@H](C)\\C=C\\CC(C)C '
               'NAME: eurysterol B sulfonic acid REASON: MISSED No steroid '
               'backbone found\n'
               ' * SMILES: '
               'C[C@]12CC[C@@H]3[C@H]([C@H]1CCC2=O)CCC4=C3C=CC(=C4)OS(=O)(=O)O '
               'NAME: sulfuric acid '
               '[(8R,9R,13S,14R)-13-methyl-17-oxo-7,8,9,11,12,14,15,16-octahydro-6H-cyclopenta[a]phenanthren-3-yl] '
               'ester REASON: MISSED No steroid backbone found\n'
               ' * SMILES: '
               'S(O[C@@H]1[C@@]2([C@]([C@]3([C@](CC2)(C4=C(CC3)C=C(OS(O)(=O)=O)C=C4)[H])[H])(CC1)[H])C)(O)(=O)=O '
               'NAME: 17-beta-Estradiol-3,17-beta-sulfate REASON: MISSED No '
               'steroid backbone found\n'
               ' * SMILES: '
               'O(S(O)(=O)=O)[C@@H]1CC=2[C@]([C@]3(CC[C@]4([C@]([C@@]3(CC2)[H])(CC[C@@]4([C@H](C)CCCC(C)C)[H])[H])C)[H])(C)CC1 '
               'NAME: cholesterol sulfate REASON: MISSED No steroid backbone '
               'found\n'
               ' * SMILES: '
               'S(O[C@@H]1[C@@]2(C(C3C([C@@]4(C(CC3)=CC(=O)CC4)C)CC2)CC1)C)(O)(=O)=O '
               'NAME: 4-Androsten-17.beta.-ol-3-one sulfate REASON: MISSED No '
               'steroid backbone found\n'
               ' * SMILES: '
               'C1=C(OS(O)(=O)=O)C=CC2=C1CC[C@]3([C@@]4(CC[C@]([C@]4(CC[C@@]32[H])C)(C#C)O)[H])[H] '
               'NAME: 17alpha-ethynylestradiol 3-sulfate REASON: MISSED No '
               'steroid backbone found\n'
               ' * SMILES: '
               'C1C2=CC(CC[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C1)[H])(CC[C@H]4OS(=O)([O-])=O)[H])C)[H])C)=O '
               'NAME: epitestosterone sulfate(1-) REASON: MISSED No steroid '
               'backbone found\n'
               ' * SMILES: '
               '[C@@H]1(OS(O)(=O)=O)[C@]2(CC[C@]3([C@]([C@@]2(CC1)[H])(CC[C@@]4([C@@]3(CC[C@@H](C4)O)C)[H])[H])[H])C '
               'NAME: (3beta,5alpha,17beta)-3-hydroxyandrostan-17-yl sulfate '
               'REASON: MISSED No steroid backbone found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O[C@H]1[C@@H]([C@@H](C(=O)C1)/C=C/[C@H](O)C/C=C\\C/C=C/C/C=C\\CC)CCC(O)=O',
                                     'name': 'ent-11-D4t-NeuroP',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C1[C@@H]2[C@@H]([C@@H](N2CC3=CC=CC=C3Cl)CO)C4=CC=CC=C4N1S(=O)(=O)C5=CC=C(C=C5)F',
                                     'name': '[(1R,2aS,8bS)-2-[(2-chlorophenyl)methyl]-4-(4-fluorophenyl)sulfonyl-1,2a,3,8b-tetrahydroazeto[2,3-c]quinolin-1-yl]methanol',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C=CC1C2C[C@]3([C@@]45C(C6=CC=CC=C6N4)[C@@](CN3)(O[C@@](C25C(=O)OC)(O[C@]1(O[C@@H]7[C@H]([C@@H]([C@H]([C@H](CO)O7)O)O)O)[H])[H])[H])[H]',
                                     'name': 'Cymoside',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'CC1(CC2=C(C(NC3=C2C4=CC=CC=C4C=C3)C5=CC(=C(C=C5)OCC(=O)O)OC)C(=O)C1)C',
                                     'name': '2-[4-(2,2-dimethyl-4-oxo-1,3,5,6-tetrahydrobenzo[a]phenanthridin-5-yl)-2-methoxyphenoxy]acetic '
                                             'acid',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'CC1(C)C2CC(=O)C1(C)CC2O',
                                     'name': '5-hydroxycamphor',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'O1C(C(=O)N2C(CCC2)C(=O)NC(C(CC)C)C(=O)N(C(C(C)C)C(=O)N(C(C(=O)NCCC1=O)C)C)C)CC(O)(C)C',
                                     'name': 'Hydroxydestruxin B',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'CCC1=NN(C(=O)C2=CC3=C(N21)C=CO3)C(C)C(=O)NC(C)CCC4=CC=CO4',
                                     'name': '2-(1-ethyl-4-oxo-3-furo[3,4]pyrrolo[3,5-c][1,2,4]triazinyl)-N-[4-(2-furanyl)butan-2-yl]propanamide',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'O[C@@H]1C(=C)[C@]2(C(C)(C)[C@H](C1)O)CC=C(CO)CC2',
                                     'name': 'Acaciicolinol D',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C[N+]1(CCCC(C1)NC(=O)CCCC(O)=O)Cc1ccc(cc1)C(=O)NCCO',
                                     'name': '3-(4-carboxybutanamido)-1-{4-[(2-hydroxyethyl)carbamoyl]benzyl}-1-methylpiperidinium',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'COC1=CC=CC=C1CNC(=O)C[C@@H]2C[C@@H]3[C@H]([C@@H](O2)CO)OC4=C3C=C(C=C4)NS(=O)(=O)C',
                                     'name': '2-[(1S,3S,4aS,9aR)-1-(hydroxymethyl)-6-(methanesulfonamido)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-[(2-methoxyphenyl)methyl]acetamide',
                                     'reason': 'No steroid backbone found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]1(CC[C@]2([H])[C@]1(C)CC[C@]1([H])[C@@]34CC[C@@H](C[C@]3(O)[C@H](O)C[C@@]21OC4)OS(O)(=O)=O)[C@H](C)\\C=C\\CC(C)C',
                                      'name': 'eurysterol B sulfonic acid',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(=O)COS(O)(=O)=O',
                                      'name': 'cortisol 21-sulfate',
                                      'reason': 'Sulfate group not attached to '
                                                'steroid backbone via ester '
                                                'linkage'},
                                  {   'smiles': 'C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)[C@H](CC(=O)[C@]4(CO)[C@H]3CC[C@]12C)OS(O)(=O)=O)[C@H]1CC(C)=C(CO)C(=O)O1',
                                      'name': 'withalongolide G',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[H][C@@]12C[C@H](O)CC[C@]1(C)C1C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CC[C@@H](O)C(CO)COS(O)(=O)=O',
                                      'name': '5beta-scymnol sulfate',
                                      'reason': 'Sulfate group not attached to '
                                                'steroid backbone via ester '
                                                'linkage'},
                                  {   'smiles': '[H][C@@]1(CC[C@]2([H])[C@]1(C)CC[C@]1([H])[C@@]34CC[C@@H](C[C@]3(O)[C@H](O)C[C@@]21OC4)OS(O)(=O)=O)[C@H](C)CCCC(C)C',
                                      'name': 'eurysterol A sulfonic acid',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'C[C@@H]([C@H]1CC[C@H]2[C@@H]3C[C@H]4O[C@]44[C@@H](O)[C@H](CC(=O)[C@]4(C)[C@H]3CC[C@]12C)OS(O)(=O)=O)[C@H]1CC(C)=C(CO)C(=O)O1',
                                      'name': '2,3-dihydro-3beta-O-sulfate '
                                              'withaferin A',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'OS(=O)(=O)OC1CCC2C(CCC3C4CCCC4CCC23)C1',
                                      'name': 'Steryl sulfate',
                                      'reason': 'Sulfate group not attached to '
                                                'steroid backbone via ester '
                                                'linkage'},
                                  {   'smiles': 'S(OCC(=O)[C@@H]1[C@@]2([C@]([C@]3([C@@]([C@@]4(C(=CC3)C[C@@H](O)CC4)C)(CC2)[H])[H])(CC1)[H])C)(O)(=O)=O',
                                      'name': '21-hydroxypregnenolone '
                                              'monosulfate',
                                      'reason': 'Sulfate group not attached to '
                                                'steroid backbone via ester '
                                                'linkage'},
                                  {   'smiles': '[C@@]12([C@]3([C@](CC[C@@]1([C@@]4(C(C[C@H](CC4)O)=CC2)C)[H])([C@](CC3)([C@H](C)CC[C@@H](C(C)C)OS(O)(=O)=O)[H])C)[H])[H]',
                                      'name': '(24S)-hydroxycholesterol '
                                              '24-sulfate',
                                      'reason': 'Sulfate group not attached to '
                                                'steroid backbone via ester '
                                                'linkage'},
                                  {   'smiles': '[H][C@@]12C[C@H](O)CC[C@]1(C)[C@@]1([H])C[C@H](O)[C@]3(C)[C@]([H])(CC[C@@]3([H])[C@]1([H])[C@H](O)C2)[C@H](C)CCCC(CO)COS(O)(=O)=O',
                                      'name': '5beta-cyprinol sulfate',
                                      'reason': 'Sulfate group not attached to '
                                                'steroid backbone via ester '
                                                'linkage'}],
    'attempt': 1,
    'reasoning': 'The previous program failed because the SMARTS pattern used '
                 'to identify the steroid backbone was too specific or '
                 'incorrect. It did not match any of the steroid sulfate '
                 'molecules provided, resulting in "No steroid backbone found" '
                 'for all cases. To improve, we can instead identify the '
                 'steroid backbone by searching for a fused ring system '
                 'consisting of four rings with sizes 6,6,6,5 and a total of '
                 '17 carbon atoms, which is characteristic of the '
                 'cyclopentanoperhydrophenanthrene nucleus found in steroids.\n'
                 '\n'
                 'We can also modify the approach to find sulfate groups '
                 'attached via an oxygen atom to the steroid backbone without '
                 'relying on a specific SMARTS pattern for the steroid '
                 'nucleus. By analyzing the ring structures and connections, '
                 'we can more reliably identify steroid sulfates.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 53,
    'num_false_positives': 68,
    'num_true_negatives': 142169,
    'num_false_negatives': 10,
    'num_negatives': None,
    'precision': 0.4380165289256198,
    'recall': 0.8412698412698413,
    'f1': 0.5760869565217391,
    'accuracy': 0.999451862262825,
    'negative_predictive_value': 0.9999296661250958}