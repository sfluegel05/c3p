"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is derived from a triterpene, having a C30 skeleton which may be rearranged
    or missing some methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a triterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count number of carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons < 24:
        return False, f"Too few carbon atoms ({num_carbons}) to be a triterpenoid"

    # Count number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 3:
        return False, f"Too few rings ({num_rings}) to be a triterpenoid"

    # Analyze fused ring systems
    # Get list of rings (as sets of atom indices)
    ring_atom_indices = ring_info.AtomRings()
    rings = [set(ring) for ring in ring_atom_indices]

    # Build adjacency list for rings
    adjacency = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if rings[i] & rings[j]:
                adjacency[i].add(j)
                adjacency[j].add(i)

    # Find fused ring systems using BFS
    fused_ring_systems = []
    visited = set()
    for i in range(len(rings)):
        if i not in visited:
            queue = [i]
            visited.add(i)
            fused_system = set(rings[i])
            while queue:
                current = queue.pop(0)
                for neighbor in adjacency[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
                        fused_system.update(rings[neighbor])
            fused_ring_systems.append(fused_system)

    # Find the largest fused ring system
    max_fused_ring_size = max(len(system) for system in fused_ring_systems) if fused_ring_systems else 0
    if max_fused_ring_size < 10:
        return False, f"Largest fused ring system contains too few atoms ({max_fused_ring_size}) to be a triterpenoid"

    return True, "Molecule meets the criteria of a triterpenoid"


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:36615',
                          'name': 'triterpenoid',
                          'definition': 'Any terpenoid derived from a '
                                        'triterpene. The term includes '
                                        'compounds in which the C30 skeleton '
                                        'of the parent triterpene has been '
                                        'rearranged or modified by the removal '
                                        'of one or more skeletal atoms '
                                        '(generally methyl groups).',
                          'parents': ['CHEBI:26873'],
                          'xrefs': ['KEGG:C06085'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               'Attempt failed: F1 score of 0.0909090909090909 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: SMILES: '
               'O1C23OC(CC(O)(C3C4C(C5(C(C6(C(CC5)C(C(OC7OCC(O)C(OC8OC(C(O)C(O)C8O)CO)C7OC9OC(C(O)C(O)C9O)CO)CC6)(C)C)C)CC4)C)(C2)C1)C)C=C(C)C '
               'NAME: Hoduloside IV REASON: CORRECT Molecule meets the '
               'criteria of a triterpenoid\n'
               ' * SMILES: '
               'O1C23C(C4(C(C5(C(CC4)C(C(O)CC5)(C)C)C)CC2)C)(CC(O)C6(C3CC(CC6OC(=O)C)(C)C)C1)C '
               'NAME: 22-Acetylpriverogenin B REASON: CORRECT Molecule meets '
               'the criteria of a triterpenoid\n'
               ' * SMILES: '
               'O1[C@@]23C([C@]4(C([C@@]5(C(CC4)C([C@@H](O)CC5)(C)C)C)C=C2)C)(CC[C@@]6(C3C([C@@H](CC6)C)C)C1=O)C '
               'NAME: DEHYDRO (11,12)URSOLIC ACID LACTONE REASON: CORRECT '
               'Molecule meets the criteria of a triterpenoid\n'
               ' * SMILES: '
               'O1C23C(C4(C(C5(C(CC4)C(C(OC(=O)C)CC5)(C)C)C)CC2O)C)(CCC6(C3CC(CC6)(C)C)C1=O)C '
               'NAME: Oleanderolide 3-acetate REASON: CORRECT Molecule meets '
               'the criteria of a triterpenoid\n'
               ' * SMILES: '
               'O1[C@]2([C@@H]([C@@]3([C@H](C(CCC3)(C)C)CC2)C)CC[C@H]4[C@]1(CC[C@@H]5[C@@]4(CCCC5(C)C)C)C)C '
               'NAME: Onoceranoxide REASON: CORRECT Molecule meets the '
               'criteria of a triterpenoid\n'
               'False positives: NONE\n'
               'False negatives: SMILES: '
               'OC1(C(CC(=O)CC1C)(C)C)CC(=O)/C(/C)=C/C=C\\C(\\C)=C\\C=C\\C=C(/C=C\\C=C(/C#CC=2C(CC(O)CC2C)(C)C)\\COC(=O)CCCCC)\\C '
               "NAME: 19'-Hexanoyloxymytiloxanthin REASON: MISSED Too few "
               'rings (2) to be a triterpenoid\n'
               ' * SMILES: '
               'O1[C@@]2([C@@]([C@@]3([C@]([C@]4([C@@]([C@@]5([C@@]([C@@H](O)C4)(C[C@@H](O)CC5)[H])C)(CC3)[H])[H])(C2)[H])C)([C@@H]([C@]16OC[C@@H](CC6)C)C)[H])[H] '
               'NAME: Chlorogenin REASON: MISSED Largest fused ring system '
               'contains too few atoms (14) to be a triterpenoid\n'
               ' * SMILES: '
               'OC1C(C2C(C3C(C4(C(C5C(CC4)(C=CC(C5)(C)C)C)=CC3)C)(CC2)C)(CC1)C)(C)C '
               'NAME: Anhydrosophoradiol REASON: MISSED Largest fused ring '
               'system contains too few atoms (14) to be a triterpenoid\n'
               ' * SMILES: '
               'CC1(C)CC[C@@]2(CC[C@]3(C)C(=CC[C@@H]4[C@@]5(C)CC[C@H](O)[C@@](C)(CO)[C@@H]5[C@H](O)C[C@@]34C)[C@@H]2C1)C(=O)O[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O '
               'NAME: 3beta,6beta,23-trihydroxyolean-12-en-28-oic acid '
               '28-O-beta-D-glucopyranoside REASON: MISSED Largest fused ring '
               'system contains too few atoms (10) to be a triterpenoid\n'
               ' * SMILES: '
               'OC1C2C(C3C(C4(C(C5C(CC4)(CCC(C5)(C)C)C(O)=O)=CC3)C)(C1)C)(CCC(=O)C2(C)C)C '
               'NAME: 6beta-Hydroxy-3-oxo-12-oleanen-28-oic acid REASON: '
               'MISSED Largest fused ring system contains too few atoms (14) '
               'to be a triterpenoid\n'
               ' * SMILES: '
               'COC(=O)C[C@H]1[C@@]2(C)C[C@@]3(OC(C)=O)[C@]1(C)[C@H]1CC[C@@]4(C)[C@@H](OC(=O)[C@H](O)C4=C1[C@H](OC(C)=O)[C@@]3(OC(C)=O)[C@H]2OC(=O)C(\\C)=C\\C)c1ccoc1 '
               'NAME: trichagmalin E REASON: MISSED Largest fused ring system '
               'contains too few atoms (10) to be a triterpenoid\n'
               ' * SMILES: '
               'O=C1C([C@H]2[C@](C=3C([C@]4([C@]([C@@H]([C@H](C(=O)O)CCCC(O)(C)C)[C@@H](C4)O)(C)[C@@H](C3)O)C)=CC2)(C)CC1)(C)C '
               'NAME: Pinicolic acid F REASON: MISSED Largest fused ring '
               'system contains too few atoms (10) to be a triterpenoid\n'
               ' * SMILES: '
               'O[C@H]1[C@@]2([C@@]([C@@]3([C@]([C@]4(C([C@]5([C@@](CC4)(CCC(C5)(C)C)C(O[C@@H]6OC[C@H](O)[C@H](O)[C@H]6O)=O)[H])=CC3)C)(C1)C)[H])(C[C@H](O)[C@H](O)[C@]2(CO)C)C)[H] '
               'NAME: Madlongiside C REASON: MISSED Largest fused ring system '
               'contains too few atoms (14) to be a triterpenoid\n'
               ' * SMILES: '
               'O1[C@@H](O[C@H]2[C@]3([C@H](C4=CCC5[C@@]6(C([C@]([C@@H](O)CC6)(CO)C)CC[C@]5([C@]4(C)CC3)C)C)CC([C@H]2O)(C)C)C)[C@H](O)[C@@H](O)[C@H](C1)O '
               'NAME: MK800-62F1 REASON: MISSED Largest fused ring system '
               'contains too few atoms (14) to be a triterpenoid\n'
               ' * SMILES: '
               'O(C1C2C3(C(C(CC3)C(CC(O)C=C(C)C)C)(CCC2(C4C(=C1)C(C(O)CC4)(C)C)C)C)C)C5OC(C(O)C(O)C5O)CO '
               'NAME: (3b,7b,22x)-Cucurbita-5,24-diene-3,7,23-triol '
               '7-glucoside REASON: MISSED Largest fused ring system contains '
               'too few atoms (10) to be a triterpenoid\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'OC1(C(CC(=O)CC1C)(C)C)CC(=O)/C(/C)=C/C=C\\C(\\C)=C\\C=C\\C=C(/C=C\\C=C(/C#CC=2C(CC(O)CC2C)(C)C)\\COC(=O)CCCCC)\\C',
                                      'name': "19'-Hexanoyloxymytiloxanthin",
                                      'reason': 'Too few rings (2) to be a '
                                                'triterpenoid'},
                                  {   'smiles': 'O1[C@@H](OC(C(O)(CC/C=C(/CCC=C(C)C)\\C)C)CCC(O)C(O)(CC/C=C(/CCC=C(C)C)\\C)C)[C@@H](O)[C@@H](O)[C@@H]([C@H]1CO)O',
                                      'name': '14,15-dihydroxyirpexan',
                                      'reason': 'Too few rings (1) to be a '
                                                'triterpenoid'},
                                  {   'smiles': 'CC(C)=CCC\\C(C)=C\\C=C\\C(C)=C\\C=C\\C=C(C)\\C=C\\C=C(/C)CCC=C(C)C',
                                      'name': "all-trans-4,4'-diapo-zeta-carotene",
                                      'reason': 'Too few rings (0) to be a '
                                                'triterpenoid'},
                                  {   'smiles': 'O1C2C3(C(C4(C(C2)C(C=C(OC)C4OC)C)C)C(=O)C(OC)=C(C3CC1=O)C)C',
                                      'name': '1alpha-O-Methylquassin',
                                      'reason': 'Too few carbon atoms (23) to '
                                                'be a triterpenoid'},
                                  {   'smiles': 'C\\C(CCCC(C)(C)O)=C/C=C/C(C)=C/C=C/C(C)=C/C=C/C=C(C)/C=C/C=C(C)/C=C/c1c(C)ccc(C)c1C',
                                      'name': 'Hydroxychlorobactene',
                                      'reason': 'Too few rings (1) to be a '
                                                'triterpenoid'},
                                  {   'smiles': 'COC1=C(C)[C@]2(O)CC(=O)O[C@@H]3C[C@H]4[C@H](C)C[C@H](O)C(=O)[C@]4(C)[C@@H](C1=O)[C@]23C',
                                      'name': 'Picrasin G',
                                      'reason': 'Too few carbon atoms (21) to '
                                                'be a triterpenoid'},
                                  {   'smiles': 'CC(=C)[C@@H]1CCC(=C)[C@H](CC[C@H]2C(C)=CC[C@@H](C(C)=C)[C@]2(C)CCC(O)=O)[C@@]1(C)CCC(O)=O',
                                      'name': 'lansic acid',
                                      'reason': 'Too few rings (2) to be a '
                                                'triterpenoid'},
                                  {   'smiles': 'O(C(=O)/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(\\C=C\\C=1C(CCCC1C)(C)C)/C)\\C)/C)/C)CC',
                                      'name': 'all-trans-Carophyll yellow',
                                      'reason': 'Too few rings (1) to be a '
                                                'triterpenoid'},
                                  {   'smiles': 'O=C(O[C@@H]1OC[C@@H](O)[C@@H]([C@H]1OC(=O)CCCCCCCCCCC(C)C)O)/C(=C/C=C/C(=C/C=C/C(=C/C=C/C=C(/C=C/C=C(/C=C/C=C(/C(=O)O)\\C)\\C)\\C)/C)/C)/C',
                                      'name': 'Diapolycopenedioic acid xylosyl '
                                              'ester',
                                      'reason': 'Too few rings (1) to be a '
                                                'triterpenoid'},
                                  {   'smiles': 'OC1CC(C(=C(C1)C)\\C=C\\C(\\C)=C/C=C/C(/C)=C\\C=C\\C=C(\\C=C/C=C(\\C)/C=C\\C(=O)C)/C)(C)C',
                                      'name': 'Reticulataxanthin',
                                      'reason': 'Too few rings (1) to be a '
                                                'triterpenoid'}],
    'attempt': 1,
    'reasoning': 'The previous program failed to correctly classify many '
                 'triterpenoids because the method for identifying fused ring '
                 'systems was insufficient. Specifically, the code attempted '
                 'to build fused ring systems by simply merging rings that '
                 'share atoms, but it did not account for all connected rings '
                 'properly. This led to underestimating the size of the '
                 'largest fused ring system in some molecules.\n'
                 '\n'
                 'Moreover, the criteria for classification may have been too '
                 'strict, particularly the requirement that the largest fused '
                 'ring system must contain at least 15 atoms. Many genuine '
                 'triterpenoids have complex ring systems that may not meet '
                 'this exact criterion but are still triterpenoids.\n'
                 '\n'
                 'To improve the classification:\n'
                 '\n'
                 '1. **Enhance fused ring detection**: Use a more robust '
                 'method to identify fused ring systems by constructing a ring '
                 'adjacency graph and finding connected components. This '
                 'ensures all connected rings are considered when determining '
                 'the size of fused ring systems.\n'
                 '\n'
                 '2. **Adjust thresholds**: Lower the minimum number of '
                 'carbons and the size of the largest fused ring system '
                 'required. This makes the criteria more inclusive, capturing '
                 'more true positives without significantly increasing false '
                 'positives.\n'
                 '\n'
                 "Here's the improved code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 1352,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.931129476584022,
    'f1': 0.964336661911555,
    'accuracy': 0.931129476584022,
    'negative_predictive_value': 0.0}