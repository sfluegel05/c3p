"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: cucurbitacin
"""

from rdkit import Chem
from rdkit.Chem import Descriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are tetracyclic triterpenoids derived from cucurbitane.
    They have a characteristic fused ring system of three six-membered rings
    and one five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cucurbitacin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # List of tuples of atom indices
    num_rings = len(atom_rings)
    if num_rings < 4:
        return False, f"Contains {num_rings} rings, less than 4 required for tetracyclic system"

    # Build a list of ring atoms and ring sizes
    ring_atoms_list = [set(ring) for ring in atom_rings]
    ring_sizes = [len(ring) for ring in atom_rings]

    # Build ring fusion graph
    # Nodes are ring indices, edges exist if rings share two or more atoms
    ring_graph = {}
    for i in range(num_rings):
        ring_graph[i] = set()
        for j in range(num_rings):
            if i != j:
                # Check if rings i and j are fused
                shared_atoms = ring_atoms_list[i] & ring_atoms_list[j]
                if len(shared_atoms) >= 2:
                    ring_graph[i].add(j)

    # Find connected components (fused ring systems)
    visited = set()
    fused_ring_systems = []

    for i in range(num_rings):
        if i not in visited:
            # Perform DFS to find all rings in this fused system
            stack = [i]
            component = set()
            while stack:
                ring_idx = stack.pop()
                if ring_idx not in visited:
                    visited.add(ring_idx)
                    component.add(ring_idx)
                    stack.extend(ring_graph[ring_idx] - visited)
            fused_ring_systems.append(component)

    # Look for a fused ring system with four rings: 3 six-membered and 1 five-membered
    tetracyclic_found = False
    for system in fused_ring_systems:
        if len(system) == 4:
            sizes = [ring_sizes[idx] for idx in system]
            num_6 = sizes.count(6)
            num_5 = sizes.count(5)
            if num_6 == 3 and num_5 == 1:
                tetracyclic_found = True
                break

    if not tetracyclic_found:
        return False, "Does not contain the characteristic tetracyclic ring system of cucurbitacins"

    # Check for triterpenoid skeleton (30 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 30:
        return False, f"Carbon count {c_count} is less than 30, not a triterpenoid"

    # Check for characteristic functional groups
    # Cucurbitacins often have multiple hydroxyl and ketone groups
    num_hydroxyl = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[OX2H]')))
    num_ketone = len(mol.GetSubstructMatches(Chem.MolFromSmarts('[$([CX3]=O)]')))
    if num_hydroxyl + num_ketone < 3:
        return False, f"Contains {num_hydroxyl} hydroxyl and {num_ketone} ketone groups, less than 3 total"

    # Passed all checks
    return True, "Contains characteristic tetracyclic ring system and functional groups of cucurbitacins"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'cucurbitacin',
        'definition': 'Any one of a class of tetracyclic triterpenoids, formally derived from the triterpene hydrocarbon cucurbitane, developed by some plants (especially those of the family Cucurbitaceaeas) as a defence mechanism against herbivores.',
        'parents': []
    }
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:16219',
                          'name': 'cucurbitacin',
                          'definition': 'Any one of a class of tetracyclic '
                                        'triterpenoids, formally derived from '
                                        'the triterpene hydrocarbon '
                                        'cucurbitane, developed by some plants '
                                        '(especially those of the family '
                                        'Cucurbitaceaeas) as a defence '
                                        'mechanism against herbivores.',
                          'parents': ['CHEBI:26893'],
                          'xrefs': ['Wikipedia:Cucurbitacin'],
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
               'O1C23C(C4(C(C5(C(C(CC5)C(C/C=C/C(OC6OC(C(O)C(O)C6O)CO)(C)C)C)(CC4)C)C)C=C2)C1OC)CCC(OC7OC(C(O)C(O)C7O)CO)C3(C)C '
               'NAME: Goyaglycoside g REASON: MISSED Does not contain the '
               'cucurbitane tetracyclic skeleton\n'
               ' * SMILES: '
               'O([C@H]1[C@H]([C@]2([C@@]([C@@]3([C@](CC2)([C@@]4(C(=C[C@@H]3O)C([C@H](O[C@@H]5OC[C@@H](OC(=O)C)[C@H](OC(=O)C)C5O)CC4)(C)C)[H])C)[H])(C1)C)C)C(CCC=C(C)C)C)[C@@H]6O[C@H]([C@@H](OC(=O)C)[C@H](O)[C@@H]6O)COC(=O)C '
               'NAME: Hebevinoside XIII REASON: MISSED Does not contain the '
               'cucurbitane tetracyclic skeleton\n'
               ' * SMILES: '
               'CC(C)(O)[C@H](O)CC(=O)[C@](C)(O)[C@H]1[C@H](O)C[C@@]2(C)[C@@H]3CC=C4[C@@H](C=C(O)C(=O)C4(C)C)[C@]3(C)C(=O)C[C@]12C '
               'NAME: cucurbitacin J REASON: MISSED Does not contain the '
               'cucurbitane tetracyclic skeleton\n'
               ' * SMILES: '
               'O1C23C(C4(C(C5(C(C(CC5)C(C/C=C/C(O)(C)C)C)(CC4)C)C)C=C2)C1OC)CCC(OC6OC(C(O)C(O)C6O)CO)C3(C)C '
               'NAME: Goyaglycoside a REASON: MISSED Does not contain the '
               'cucurbitane tetracyclic skeleton\n'
               ' * SMILES: '
               'O([C@@H]1[C@]2([C@]3([C@@]([C@](CC3)([C@@H](C/C=C/C(O)(C)C)C)[H])(CC[C@]2([C@]4(C(=C1)C([C@@H](O)CC4)(C)C)[H])C=O)C)C)[H])[C@@H]5O[C@@H]([C@@H](O)[C@H](O)[C@H]5O)CO '
               'NAME: Momordicoside L REASON: MISSED Does not contain the '
               'cucurbitane tetracyclic skeleton\n'
               ' * SMILES: '
               'O(C1C(C2(C(C3C(CC2)(C4C(=CC3)C(C(O)C(OC5OC(C(O)C(O)C5O)CO)C4)(C)C)C)(C1)C)C)C(O)(C)C(=O)CCC(O)(C)C)C(=O)C '
               'NAME: '
               '[17-(2,6-Dihydroxy-6-methyl-3-oxoheptan-2-yl)-3-hydroxy-4,4,9,13,14-pentamethyl-2-[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxy-2,3,7,8,10,11,12,15,16,17-decahydro-1H-cyclopenta[a]phenanthren-16-yl] '
               'acetate REASON: MISSED Does not contain the cucurbitane '
               'tetracyclic skeleton\n'
               ' * SMILES: '
               'CC(C)(O)CCC(=O)[C@](C)(O)[C@H]1[C@H](O)C[C@@]2(C)[C@@H]3CC=C4[C@@H](C[C@@H](O)[C@@H](O)C4(C)C)[C@]3(C)C(=O)C[C@]12C '
               'NAME: Cucurbitacin P REASON: MISSED Does not contain the '
               'cucurbitane tetracyclic skeleton\n'
               ' * SMILES: '
               'C[C@H](CC\\C=C(/C)CO[C@@H]1O[C@H](CO[C@@H]2O[C@H](CO)[C@@H](O)[C@H](O)[C@H]2O)[C@@H](O)[C@H](O)[C@H]1O)[C@H]1CC[C@@]2(C)[C@@H]3CC=C4[C@@H](CC[C@H](O)C4(C)C)[C@]3(C)C(=O)C[C@]12C '
               'NAME: Carnosifloside I REASON: MISSED Does not contain the '
               'cucurbitane tetracyclic skeleton\n'
               ' * SMILES: '
               'O1[C@@]23O[C@@]([C@]1(CC(C3[C@@]4([C@](C5[C@]6([C@]7(C6)C(CC5)C([C@@H](O[C@H]8OC[C@H](O)[C@@H](O)[C@@H]8O)CC7)(C)C)CC4)(C2O)C)C)C)[H])(C(O)(C)C)[H] '
               'NAME: '
               '(2R,3S,4R,5S)-2-[[(1S,3S,9S,12R,14S,17R,19R,21R,22S)-2-hydroxy-22-(2-hydroxypropan-2-yl)-3,8,8,17,19-pentamethyl-23,24-dioxaheptacyclo[19.2.1.01,18.03,17.04,14.07,12.012,14]tetracosan-9-yl]oxy]oxane-3,4,5-triol '
               'REASON: MISSED Does not contain the cucurbitane tetracyclic '
               'skeleton\n'
               ' * SMILES: '
               'O([C@@H]1[C@@]([C@@]2([C@](C3[C@@](CC2)([C@@]4(C(=C[C@H]3O)C(C(O[C@@H]5OC[C@@H](OC(=O)C)[C@H](O)[C@H]5O)CC4)(C)C)[H])C)(C1)C)C)([C@@H](CCC=C(C)C)C)[H])[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)COC(=O)C '
               'NAME: Hebevinoside II REASON: MISSED Does not contain the '
               'cucurbitane tetracyclic skeleton\n'
               '------\n'
               '\n'
               'In your reasoning step, analyze the previous program and the '
               'above outcomes, hypothesizing about what went wrong, and how '
               'to improve.\n',
    'sample_true_negatives': [   {   'smiles': 'O[C@H]1[C@@H]([C@@H](C(=O)C1)/C=C/[C@H](O)C/C=C\\C/C=C/C/C=C\\CC)CCC(O)=O',
                                     'name': 'ent-11-D4t-NeuroP',
                                     'reason': 'Contains 1 rings, less than 4 '
                                               'required for tetracyclic '
                                               'system'},
                                 {   'smiles': 'C1[C@@H]2[C@@H]([C@@H](N2CC3=CC=CC=C3Cl)CO)C4=CC=CC=C4N1S(=O)(=O)C5=CC=C(C=C5)F',
                                     'name': '[(1R,2aS,8bS)-2-[(2-chlorophenyl)methyl]-4-(4-fluorophenyl)sulfonyl-1,2a,3,8b-tetrahydroazeto[2,3-c]quinolin-1-yl]methanol',
                                     'reason': 'Does not contain the '
                                               'characteristic tetracyclic '
                                               'ring system of cucurbitacins'},
                                 {   'smiles': 'C=CC1C2C[C@]3([C@@]45C(C6=CC=CC=C6N4)[C@@](CN3)(O[C@@](C25C(=O)OC)(O[C@]1(O[C@@H]7[C@H]([C@@H]([C@H]([C@H](CO)O7)O)O)O)[H])[H])[H])[H]',
                                     'name': 'Cymoside',
                                     'reason': 'Does not contain the '
                                               'characteristic tetracyclic '
                                               'ring system of cucurbitacins'},
                                 {   'smiles': 'CC1(CC2=C(C(NC3=C2C4=CC=CC=C4C=C3)C5=CC(=C(C=C5)OCC(=O)O)OC)C(=O)C1)C',
                                     'name': '2-[4-(2,2-dimethyl-4-oxo-1,3,5,6-tetrahydrobenzo[a]phenanthridin-5-yl)-2-methoxyphenoxy]acetic '
                                             'acid',
                                     'reason': 'Does not contain the '
                                               'characteristic tetracyclic '
                                               'ring system of cucurbitacins'},
                                 {   'smiles': 'CC1(C)C2CC(=O)C1(C)CC2O',
                                     'name': '5-hydroxycamphor',
                                     'reason': 'Contains 2 rings, less than 4 '
                                               'required for tetracyclic '
                                               'system'},
                                 {   'smiles': 'O1C(C(=O)N2C(CCC2)C(=O)NC(C(CC)C)C(=O)N(C(C(C)C)C(=O)N(C(C(=O)NCCC1=O)C)C)C)CC(O)(C)C',
                                     'name': 'Hydroxydestruxin B',
                                     'reason': 'Contains 2 rings, less than 4 '
                                               'required for tetracyclic '
                                               'system'},
                                 {   'smiles': 'CCC1=NN(C(=O)C2=CC3=C(N21)C=CO3)C(C)C(=O)NC(C)CCC4=CC=CO4',
                                     'name': '2-(1-ethyl-4-oxo-3-furo[3,4]pyrrolo[3,5-c][1,2,4]triazinyl)-N-[4-(2-furanyl)butan-2-yl]propanamide',
                                     'reason': 'Does not contain the '
                                               'characteristic tetracyclic '
                                               'ring system of cucurbitacins'},
                                 {   'smiles': 'O[C@@H]1C(=C)[C@]2(C(C)(C)[C@H](C1)O)CC=C(CO)CC2',
                                     'name': 'Acaciicolinol D',
                                     'reason': 'Contains 2 rings, less than 4 '
                                               'required for tetracyclic '
                                               'system'},
                                 {   'smiles': 'C[N+]1(CCCC(C1)NC(=O)CCCC(O)=O)Cc1ccc(cc1)C(=O)NCCO',
                                     'name': '3-(4-carboxybutanamido)-1-{4-[(2-hydroxyethyl)carbamoyl]benzyl}-1-methylpiperidinium',
                                     'reason': 'Contains 2 rings, less than 4 '
                                               'required for tetracyclic '
                                               'system'},
                                 {   'smiles': 'COC1=CC=CC=C1CNC(=O)C[C@@H]2C[C@@H]3[C@H]([C@@H](O2)CO)OC4=C3C=C(C=C4)NS(=O)(=O)C',
                                     'name': '2-[(1S,3S,4aS,9aR)-1-(hydroxymethyl)-6-(methanesulfonamido)-3,4,4a,9a-tetrahydro-1H-pyrano[3,4-b]benzofuran-3-yl]-N-[(2-methoxyphenyl)methyl]acetamide',
                                     'reason': 'Does not contain the '
                                               'characteristic tetracyclic '
                                               'ring system of cucurbitacins'}],
    'sample_false_negatives': [   {   'smiles': 'O1C23C(C4(C(C5(C(C(CC5)C(C/C=C/C(OC6OC(C(O)C(O)C6O)CO)(C)C)C)(CC4)C)C)C=C2)C1OC)CCC(OC7OC(C(O)C(O)C7O)CO)C3(C)C',
                                      'name': 'Goyaglycoside g',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': 'O1C23C(C4(C(C5(C(C(CC5)C(C/C=C/C(O)(C)C)C)(CC4)C)C)C=C2)C1OC)CCC(OC6OC(C(O)C(O)C6O)CO)C3(C)C',
                                      'name': 'Goyaglycoside a',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': 'O1[C@@]23O[C@@]([C@]1(CC(C3[C@@]4([C@](C5[C@]6([C@]7(C6)C(CC5)C([C@@H](O[C@H]8OC[C@H](O)[C@@H](O)[C@@H]8O)CC7)(C)C)CC4)(C2O)C)C)C)[H])(C(O)(C)C)[H]',
                                      'name': '(2R,3S,4R,5S)-2-[[(1S,3S,9S,12R,14S,17R,19R,21R,22S)-2-hydroxy-22-(2-hydroxypropan-2-yl)-3,8,8,17,19-pentamethyl-23,24-dioxaheptacyclo[19.2.1.01,18.03,17.04,14.07,12.012,14]tetracosan-9-yl]oxy]oxane-3,4,5-triol',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': '[H][C@@]12C[C@@H](C)[C@]3([H])[C@@]4(C)CC[C@@]56C[C@@]55CC[C@H](O[C@]7([H])OC[C@H](O)[C@H](O)[C@H]7O)C(C)(C)[C@]5([H])CC[C@@]6([H])[C@]4(C)[C@@H](O)[C@]3(O[C@@H]1C(C)(C)OC(C)=O)O2',
                                      'name': '25-O-acetylcimigenol '
                                              '3-O-alpha-L-arabinopyranoside',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': 'OC1C23C4(C2)C(C5(C(CC4)(C(CC5)C(CC/C=C(\\COC6OC(C(O)C(O)C6O)CO)/C)CO)C)C)CCC3C(C(OC7OC(C(O)C(O)C7O)CO)C1)(C)C',
                                      'name': '(24Z)-26-(Hexopyranosyloxy)-1,21-dihydroxy-9,19-cyclolanost-24-en-3-yl '
                                              'hexopyranoside',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': 'O1[C@@]23[C@]([C@]4([C@]([C@]5([C@@]([C@](CC5)([C@@H](C/C=C/C(O)(C)C)C)[H])(CC4)C)C)(C=C2)[H])C1)(CC[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)C3(C)C)[H]',
                                      'name': 'Momordicoside I',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': '[H][C@@]12C[C@@]3(C)[C@]4([H])CC=C5[C@@]([H])(C=C(O[C@@H]6O[C@H](CO)[C@@H](O)[C@H](O)[C@H]6O)C(=O)C5(C)C)[C@]4(C)C(=O)C[C@]3(C)[C@@]1([H])[C@@](C)(O)C(=O)C[C@H](O2)C(C)(C)O',
                                      'name': '(16alpha,20R,24S)-16,24-epoxy-2,20,25-trihydroxycucurbita-1,5-diene-3,11,22-trione '
                                              '2-O-beta-D-glucopyranoside',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': 'O1[C@@]23[C@]([C@]4([C@]([C@]5([C@@]([C@](CC5)([C@@H](C/C=C/C(OC)(C)C)C)[H])(CC4)C)C)(C=C2)[H])C1OC)(CC[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@H](O)[C@H]6O)CO)C3(C)C)[H]',
                                      'name': 'Goyaglycoside c',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': 'O1[C@@]23[C@]([C@]4([C@]([C@]5([C@@]([C@](CC5)([C@@H](C/C=C/C(OC)(C)C)C)[H])(CC4)C)C)(C=C2)[H])C1)(CC[C@H](O[C@@H]6O[C@@H]([C@@H](O)[C@@H](O)[C@H]6O)CO)C3(C)C)[H]',
                                      'name': 'Momordicoside G',
                                      'reason': 'Does not contain the '
                                                'characteristic tetracyclic '
                                                'ring system of cucurbitacins'},
                                  {   'smiles': 'OC1C(C=2C(C3(C(C4(C(C(CC4)C(CCC=C(C)C)C)(CC3)C)C)CC2)C)CC1)(C)C',
                                      'name': '3beta-Cucurbita-5,24-dien-3-ol',
                                      'reason': 'Contains 1 hydroxyl and 0 '
                                                'ketone groups, less than 3 '
                                                'total'}],
    'attempt': 1,
    'reasoning': 'After analyzing the previous code and the outcomes, I '
                 'realized that the SMARTS pattern used to detect the '
                 "cucurbitane skeleton was too specific and didn't match the "
                 'actual structures of cucurbitacins. This led to false '
                 'negatives, as the pattern failed to identify the cucurbitane '
                 'tetracyclic core present in these molecules.\n'
                 '\n'
                 "To improve, I'll develop a method to detect the presence of "
                 'a fused tetracyclic ring system consisting of three '
                 'six-membered rings and one five-membered ring fused '
                 'together, which is characteristic of cucurbitacins. By '
                 'analyzing the ring systems and their connections, we can '
                 'more accurately identify cucurbitacins.\n'
                 '\n'
                 "Here's the updated code:",
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 57,
    'num_false_positives': 100,
    'num_true_negatives': 9696,
    'num_false_negatives': 10,
    'num_negatives': None,
    'precision': 0.3630573248407643,
    'recall': 0.8507462686567164,
    'f1': 0.5089285714285715,
    'accuracy': 0.9888472067322316,
    'negative_predictive_value': 0.9989697094580672}