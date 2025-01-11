"""
Classifies: CHEBI:26125 phytosterols
"""
"""
Classifies: phytosterols
"""
from rdkit import Chem

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    Phytosterols are sterols similar to cholesterol which occur in plants and 
    vary only in carbon side chains and/or presence or absence of a double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phytosterol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Build ring adjacency list to identify fused rings
    ring_adj = [set() for _ in rings]
    for i, ring1 in enumerate(rings):
        for j, ring2 in enumerate(rings):
            if i >= j:
                continue
            if set(ring1) & set(ring2):
                # Rings i and j are fused (share at least one atom)
                ring_adj[i].add(j)
                ring_adj[j].add(i)

    # Find fused ring systems using DFS
    visited = set()
    fused_ring_systems = []
    for i in range(len(rings)):
        if i in visited:
            continue
        stack = [i]
        fused_rings = set()
        while stack:
            ring_idx = stack.pop()
            if ring_idx in visited:
                continue
            visited.add(ring_idx)
            fused_rings.add(ring_idx)
            # Add neighboring rings that are fused
            stack.extend(ring_adj[ring_idx] - visited)
        fused_ring_systems.append(fused_rings)

    # Look for steroid backbone: a fused ring system with 4 rings,
    # consisting of 3 six-membered rings and 1 five-membered ring
    has_steroid_backbone = False
    steroid_ring_system = None
    for ring_system in fused_ring_systems:
        if len(ring_system) != 4:
            continue
        ring_sizes = [len(rings[ring_idx]) for ring_idx in ring_system]
        num_six = ring_sizes.count(6)
        num_five = ring_sizes.count(5)
        if num_six == 3 and num_five == 1:
            has_steroid_backbone = True
            steroid_ring_system = ring_system
            break

    if not has_steroid_backbone:
        return False, "No steroid backbone found"

    # Identify atoms in the steroid backbone
    backbone_atoms = set()
    for ring_idx in steroid_ring_system:
        backbone_atoms.update(rings[ring_idx])

    # Identify side chains attached to the backbone
    side_chain_atoms = set()
    atoms_to_visit = set()
    for atom_idx in backbone_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx not in backbone_atoms and nbr_idx not in side_chain_atoms:
                side_chain_atoms.add(nbr_idx)
                atoms_to_visit.add(nbr_idx)

    # Traverse side chain atoms to get full side chain
    visited_atoms = set()
    while atoms_to_visit:
        current_idx = atoms_to_visit.pop()
        if current_idx in visited_atoms:
            continue
        visited_atoms.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        for neighbor in current_atom.GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx not in backbone_atoms and nbr_idx not in visited_atoms:
                side_chain_atoms.add(nbr_idx)
                atoms_to_visit.add(nbr_idx)

    # Count the number of carbon atoms in side chains
    num_side_chain_carbons = sum(1 for idx in side_chain_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)

    # Cholesterol has 8 carbons in the side chain; phytosterols typically have more
    if num_side_chain_carbons <= 8:
        return False, f"Side chain too short ({num_side_chain_carbons} carbons), possibly cholesterol"

    return True, "Molecule is a phytosterol with steroid backbone and appropriate side chain length"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'phytosterols',
        'definition': 'Sterols similar to cholesterol which occur in plants and vary only in carbon side chains and/or presence or absence of a double bond.',
        'parents': []
    },
    'config': {
        # Configuration details can be added here if needed
    },
    'message': None,
    'attempt': 3,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    # Additional metadata fields can be added here
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:26125',
                          'name': 'phytosterols',
                          'definition': 'Sterols similar to cholesterol which '
                                        'occur in plants and vary only in '
                                        'carbon side chains and/or presence or '
                                        'absence of a double bond.',
                          'parents': ['CHEBI:15889', 'CHEBI:26124'],
                          'xrefs': ['Wikipedia:Phytosterol'],
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
               '[H][C@@]12CC[C@@]3([H])[C@]4(C)CC[C@H]([C@H](C)CCC(=C)C(C)C)[C@@]4(C)CC[C@]33C[C@]13CC[C@H](O)[C@H]2C '
               'NAME: cycloeucalenol REASON: MISSED No steroid backbone found\n'
               ' * SMILES: '
               'OC1CC2[C@@](C3C(C4[C@@](C(CC4)[C@@H](CCCC(C)C)C)(CC3)C)CC2)(CC1)C '
               'NAME: cholestan-3-ol REASON: MISSED No steroid backbone found\n'
               ' * SMILES: '
               '[C@]1([C@@H](CCC(C(C)C)CC)C)([C@]2(CC[C@@]3([C@]4(CC[C@@H](C[C@]4(CC[C@]3([C@@]2(CC1)[H])[H])[H])O)C)[H])C)[H] '
               'NAME: 24-ethylcoprostanol REASON: MISSED No steroid backbone '
               'found\n'
               ' * SMILES: '
               '[H][C@@]1(O[C@H]2CC[C@@]3(C)C(C2)=CC[C@]2([H])[C@]3([H])CC[C@]3(C)[C@]([H])(CC[C@@]23[H])[C@H](C)\\C=C\\[C@@H](CC)C(C)C)O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: stigmasterol 3-O-beta-D-glucoside REASON: MISSED No '
               'steroid backbone found\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)\\C=C\\[C@@H](CC)C(C)C '
               'NAME: stigmasterol REASON: MISSED No steroid backbone found\n'
               ' * SMILES: '
               '[C@]1([C@@H](CCCC(C)C)C)([C@]2(CC[C@@]3([C@]4(CC[C@@H](C[C@]4(CC[C@]3([C@@]2(CC1)[H])[H])[H])O)C)[H])C)[H] '
               'NAME: coprostanol REASON: MISSED No steroid backbone found\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@H](C)C(C)C '
               'NAME: ergosta-5,7-dien-3beta-ol REASON: MISSED No steroid '
               'backbone found\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@@H](C)C(C)C '
               'NAME: campesterol REASON: MISSED No steroid backbone found\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@H](C)C(C)C '
               'NAME: 24-epicampesterol REASON: MISSED No steroid backbone '
               'found\n'
               ' * SMILES: '
               '[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC\\C(=C\\C)C(C)C '
               'NAME: isofucosterol REASON: MISSED No steroid backbone found\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O=C(NC1=C(C(=O)CCNC(=O)CCC(NC(=O)C)C(=O)OCC(O)CO)C=CC=C1)C(C=C)(C)C',
                                     'name': 'Citrinamide B',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'S(CC[C@H](NC(=O)[C@@H](NC(=O)[C@@H](N)CC=1C=2C(NC1)=CC=CC2)CC(O)=O)C(O)=O)C',
                                     'name': 'Trp-Asp-Met',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C(CC)(O)(C)C',
                                     'name': '2-methylbutan-2-ol',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'O1[C@@H]([C@H](O)[C@@H](O)[C@@H](NC(=O)C)C1O)CO',
                                     'name': 'N-Acetyl-D-Gulosamine',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': '[H][C@@]12CCC3=CC(=O)CC[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])CC[C@]1(O)C(C)O',
                                     'name': '17,20-dihydroxypregn-4-en-3-one',
                                     'reason': 'Side chain too short (4 '
                                               'carbons), possibly '
                                               'cholesterol'},
                                 {   'smiles': 'O=C1N2[C@H](C(=O)NCC(=O)N[C@H](C(=O)NCC(N[C@H](C(NCC(N[C@H]1CC=3C4=C(C=CC=C4)NC3)=O)=O)C(O)C5=CC(=C(O)C=C5)CC=C(C)C)=O)C(C)C)CCC2',
                                     'name': 'WIN 66306',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'CCCCCCCCCCc1c(C)c(O)c(OC)c(OC)c1O',
                                     'name': '6-decylubiquinol',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'COC1=CC(=C(C=C1Br)OC)S(=O)(=O)N2CCCCCC2',
                                     'name': '1-(4-bromo-2,5-dimethoxyphenyl)sulfonylazepane',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'C[C@@H]1CN(C(=O)C2=C(C=CC(=C2)N(C)C)O[C@@H]1CN(C)C(=O)NC3=CC=C(C=C3)F)[C@@H](C)CO',
                                     'name': '1-[[(2S,3R)-8-(dimethylamino)-5-[(2S)-1-hydroxypropan-2-yl]-3-methyl-6-oxo-3,4-dihydro-2H-1,5-benzoxazocin-2-yl]methyl]-3-(4-fluorophenyl)-1-methylurea',
                                     'reason': 'No steroid backbone found'},
                                 {   'smiles': 'P(=O)(N=C(N(CC)CC)C)(OCC)F',
                                     'name': 'A-234 nerve agent',
                                     'reason': 'No steroid backbone found'}],
    'sample_false_negatives': [   {   'smiles': '[H][C@@]12CC[C@@]3([H])[C@]4(C)CC[C@H]([C@H](C)CCC(=C)C(C)C)[C@@]4(C)CC[C@]33C[C@]13CC[C@H](O)[C@H]2C',
                                      'name': 'cycloeucalenol',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': 'C1[C@]2([C@]3(C[C@H](C1)O)[C@H]([C@H](C4=C2CC[C@]5([C@@]4(C(C[C@@H]5[C@@H](/C=C/[C@@H](C(C)C)C)C)=O)O)C)O)O3)C',
                                      'name': '(22E,24R)-5alpha,6alpha-epoxy3beta,7alpha,14beta-trihydroxy-ergosta-8,22,dien-15-one',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[C@]123[C@@]4([C@]([C@@]([C@@H](O)CC4)(C=O)C)(CC[C@]1([C@]5(C)CC[C@@]([C@@]5(C)CC2)([C@@H](CCC(C(C)C)=C)C)[H])[H])[H])C3',
                                      'name': '3beta-hydroxy-24-methylene-9beta-9,19-cyclolanostan-28-al',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[C@]123[C@@]4([C@](C([C@@H](O)CC4)(C)C)(CC[C@]1([C@]5([C@@]([C@@]([C@@H](CCC=C(C)C)C)(CC5)[H])(CC2)C)C)[H])[H])C3',
                                      'name': 'cycloartenol',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[C@]123[C@@]4([C@]([C@]([C@@H](O)CC4)(C)C(=O)O)(CC[C@]1([C@]5(C)CC[C@@]([C@]5(CC2)C)([C@@H](CCC=C(C)C)C)[H])[H])[H])C3',
                                      'name': '3beta-hydroxy-9beta-9,19-cyclolanost-24-en-28-oic '
                                              'acid',
                                      'reason': 'No steroid backbone found'},
                                  {   'smiles': '[C@]123[C@]([C@]4([C@]([C@](CC4)([C@@H](CCC=C(C)C)C)[H])(C)CC1)C)(CC[C@@]5([C@@]2(CC[C@@H]([C@]5(CO)C)O)C3)[H])[H]',
                                      'name': '9beta-9,19-cyclolanost-24-ene-3beta,28-diol',
                                      'reason': 'No steroid backbone found'}],
    'attempt': 3,
    'reasoning': 'The previous code failed to detect the steroid backbone in '
                 'the test molecules because the SMARTS pattern used for the '
                 'steroid backbone was too specific or incorrect. The pattern '
                 'did not account for the variations and conformations present '
                 'in phytosterols, leading to no matches. To improve, we can '
                 'implement a more general approach by analyzing the ring '
                 'systems in the molecule:\n'
                 '\n'
                 '- Identify all fused ring systems in the molecule.\n'
                 '- Look for a ring system consisting of exactly four fused '
                 'rings.\n'
                 '- Ensure that these rings include three six-membered rings '
                 'and one five-membered ring, characteristic of the steroid '
                 'backbone.\n'
                 '\n'
                 'By adopting this method, we can more accurately detect the '
                 'steroid backbone in phytosterols.',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 19,
    'num_false_positives': 100,
    'num_true_negatives': 4686,
    'num_false_negatives': 6,
    'num_negatives': None,
    'precision': 0.15966386554621848,
    'recall': 0.76,
    'f1': 0.26388888888888884,
    'accuracy': 0.9779671585948867,
    'negative_predictive_value': 0.9987212276214834}