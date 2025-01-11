"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies: sterol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol based on its SMILES string.
    A sterol is any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol
    (additional carbon atoms may be present in the side chain).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get the ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Must have at least 4 rings
    if len(rings) < 4:
        return False, f"Contains only {len(rings)} rings, requires at least 4"

    # Build ring adjacency graph
    ring_count = len(rings)
    ring_adjacency = [set() for _ in range(ring_count)]

    # Rings are connected if they share at least 2 atoms (fused)
    for i in range(ring_count):
        for j in range(i+1, ring_count):
            shared_atoms = set(rings[i]) & set(rings[j])
            if len(shared_atoms) >= 2:
                ring_adjacency[i].add(j)
                ring_adjacency[j].add(i)

    # Find fused ring systems
    visited = set()
    fused_ring_systems = []

    for i in range(ring_count):
        if i not in visited:
            stack = [i]
            component = []
            while stack:
                ring_idx = stack.pop()
                if ring_idx not in visited:
                    visited.add(ring_idx)
                    component.append(ring_idx)
                    stack.extend(ring_adjacency[ring_idx] - visited)
            fused_ring_systems.append(component)

    # Check if there is a fused ring system with at least 4 rings
    steroid_ring_system = None
    for system in fused_ring_systems:
        if len(system) >= 4:
            # Check if ring sizes are correct (three 6-membered and one 5-membered ring)
            ring_sizes = sorted([len(rings[i]) for i in system])
            if ring_sizes[:4] == [5,6,6,6]:
                steroid_ring_system = system
                break

    if steroid_ring_system is None:
        return False, "Does not contain steroid nucleus with fused 5,6,6,6 rings"

    # Define SMARTS pattern for 3-hydroxy group attached to ring A (first ring)
    hydroxy_pattern = Chem.MolFromSmarts('[#6]-1([#8H])-[#6]-[#6]-[#6]-[#6]-1')  # 3-hydroxy group
    if hydroxy_pattern is None:
        return False, "Invalid hydroxy group SMARTS pattern"

    if not mol.HasSubstructMatch(hydroxy_pattern):
        return False, "Does not have 3-hydroxy group at position 3"

    return True, "Contains steroid nucleus with 3-hydroxy group characteristic of sterols"

__metadata__ = {
    'chemical_class': {
        'name': 'sterol',
        'definition': 'Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol (additional carbon atoms may be present in the side chain).'
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:15889',
                          'name': 'sterol',
                          'definition': 'Any 3-hydroxy steroid whose skeleton '
                                        'is closely related to cholestan-3-ol '
                                        '(additional carbon atoms may be '
                                        'present in the side chain).',
                          'parents': ['CHEBI:36834'],
                          'xrefs': [   'KEGG:C00370',
                                       'LIPID_MAPS_class:LMST01',
                                       'MetaCyc:Sterols',
                                       'Wikipedia:Sterol'],
                          'all_positive_examples': []},
    'config': None,
    'message': '\n'
               "Error: module 'rdkit.Chem' has no attribute 'FusedRingUtils'\n"
               'Attempt failed: F1 score of 0 is too low.\n'
               'Outcomes:\n'
               '------\n'
               '\n'
               'True positives: NONE\n'
               'False positives: NONE\n'
               'False negatives: NONE\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': 'C1[C@]2([C@@]3(CC[C@]4([C@](C3=CC=C2C[C@H](C1)O)(CC[C@@]4([C@@H](CCCC(C)C)C)[H])[H])C)[H])C',
                                      'name': 'lumisterol 3',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCCC(C)(C)O',
                                      'name': '7alpha,25-dihydroxycholesterol',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2(C)C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](O)C(C)(C)[C@]1([H])CC3)[C@H](C)CC\\C=C(/C)CO',
                                      'name': '26-hydroxylanosterol',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': '[H][C@@]1(CC[C@@]2([H])[C@]3([H])[C@@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCCC(C)C',
                                      'name': '7beta-hydroxycholesterol',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': '[H][C@@]1(CC=C2C3=C(CC[C@]12C)[C@@]1(C)CC[C@H](O)C(C)(C)[C@]1([H])CC3)[C@H](C)CCC=C(C)C',
                                      'name': '4,4-dimethyl-5alpha-cholesta-8,14,24-trien-3beta-ol',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': 'C[C@@H](CO)CCC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3[C@H](O)C=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C',
                                      'name': '(25R)-7alpha,26-dihydroxycholesterol',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C(C=C2C[C@H](C1)O)=O)[H])(CC[C@@]4([C@H](C)CCC[C@@H](C)CO)[H])[H])C)[H])C',
                                      'name': '(25R)-3beta,26-dihydroxycholest-5-en-7-one',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': 'C1[C@@]([C@@]2([C@@](C1)(C3=CC[C@@]4([H])[C@H](C=O)[C@@H](O)CC[C@]4(C)[C@]3(CC2)[H])[H])C)([C@H](C)CCC(=C)C(C)C)[H]',
                                      'name': '4alpha-formyl-ergosta-7,24(28)-dien-3beta-ol',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': 'C1[C@@]2([C@]3(CC[C@]4([C@]([C@@]3(C(C=C2C[C@H](C1)O)=O)[H])(CC[C@@]4([C@H](C)CCC[C@@H](C)C([H])=O)[H])[H])C)[H])C',
                                      'name': '(25R)-3beta-hydroxycholest-5-en-7-one-26-al',
                                      'reason': 'Does not have 3-hydroxy group '
                                                'at position 3'},
                                  {   'smiles': '[C@]123[C@@]4([C@](C([C@@H](O)CC4)(C)C)(CC[C@]1([C@]5([C@@]([C@@]([C@@H](CCC=C(C)C)C)(CC5)[H])(CC2)C)C)[H])[H])C3',
                                      'name': 'cycloartenol',
                                      'reason': 'Does not contain steroid '
                                                'nucleus with fused 5,6,6,6 '
                                                'rings'}],
    'attempt': 2,
    'reasoning': 'The previous code failed because it attempted to use '
                 '`Chem.FusedRingUtils.IsFused(core)`, but the '
                 '`FusedRingUtils` module does not exist in the `rdkit.Chem` '
                 'package, leading to an AttributeError. This caused the '
                 'program to fail at runtime.\n'
                 '\n'
                 'Additionally, the classification had an F1 score of 0, '
                 'indicating that the code did not correctly identify any '
                 'sterols. This suggests that the method used to identify the '
                 'steroid nucleus and the 3-hydroxy group was not effective.\n'
                 '\n'
                 'To improve the code, we can:\n'
                 '\n'
                 '1. Remove the use of `Chem.FusedRingUtils` and instead '
                 'manually check if the rings are fused by analyzing shared '
                 'atoms between rings.\n'
                 '2. Use SMARTS patterns to more accurately identify the '
                 'steroid nucleus—a fused ring system consisting of three '
                 'six-membered rings and one five-membered ring.\n'
                 '3. Ensure that we are correctly identifying the 3-hydroxy '
                 'group attached to the A-ring (position 3) of the steroid '
                 'nucleus.\n'
                 '4. Account for variations in the side chain, as sterols may '
                 'have additional carbon atoms in the side chain.\n'
                 '\n'
                 'Here is the improved code:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 5,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.047619047619047616,
    'f1': 0.0909090909090909,
    'accuracy': 0.047619047619047616,
    'negative_predictive_value': 0.0}