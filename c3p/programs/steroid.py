"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: steroid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    A steroid is defined as compounds based on the cyclopenta[a]phenanthrene carbon skeleton,
    which is a tetracyclic fused ring system consisting of three six-membered rings and one five-membered ring.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    num_rings = ri.NumRings()

    # Steroids typically have 4 rings
    if num_rings < 4:
        return False, f"Contains {num_rings} rings, less than 4 rings required for steroid backbone"

    # Get ring sizes
    ring_sizes = [len(r) for r in ri.AtomRings()]

    # Check for three six-membered rings and one five-membered ring
    num_six_membered = ring_sizes.count(6)
    num_five_membered = ring_sizes.count(5)

    if num_six_membered < 3 or num_five_membered < 1:
        return False, f"Ring sizes do not match steroid backbone requirements (need at least three six-membered rings and one five-membered ring)"

    # Check if rings are fused together into one ring system
    # Create a list of sets of atoms in each ring
    ring_atom_sets = [set(r) for r in ri.AtomRings()]

    # Find fused ring systems
    fused_ring_groups = []
    rings_remaining = set(range(len(ring_atom_sets)))
    while rings_remaining:
        current_ring = rings_remaining.pop()
        fused_group = set([current_ring])
        atoms_in_group = ring_atom_sets[current_ring].copy()
        rings_to_check = set()
        for idx in rings_remaining:
            if ring_atom_sets[idx] & atoms_in_group:
                rings_to_check.add(idx)
        while rings_to_check:
            idx = rings_to_check.pop()
            rings_remaining.discard(idx)
            fused_group.add(idx)
            atoms_in_group.update(ring_atom_sets[idx])
            for idx2 in rings_remaining:
                if ring_atom_sets[idx2] & atoms_in_group:
                    rings_to_check.add(idx2)
        fused_ring_groups.append(fused_group)

    # Check if any fused ring group contains at least 4 rings
    largest_fused_group = max(fused_ring_groups, key=len)
    if len(largest_fused_group) < 4:
        return False, "No fused ring system with at least 4 rings found"

    # Check if the fused ring system contains the required ring sizes
    fused_ring_sizes = [ring_sizes[i] for i in largest_fused_group]
    if fused_ring_sizes.count(6) < 3 or fused_ring_sizes.count(5) < 1:
        return False, "Fused ring system does not have required ring sizes for steroid backbone"

    # All criteria met, molecule is classified as steroid
    return True, "Molecule contains steroid backbone of fused rings (three six-membered rings and one five-membered ring)"


__metadata__ = {  
   'chemical_class': {   'id': None,
                         'name': 'steroid',
                         'definition': 'Any of naturally occurring compounds and synthetic analogues, based on the cyclopenta[a]phenanthrene carbon skeleton, partially or completely hydrogenated; there are usually methyl groups at C-10 and C-13, and often an alkyl group at C-17. By extension, one or more bond scissions, ring expansions and/or ring contractions of the skeleton may have occurred. Natural steroids are derived biogenetically from squalene which is a triterpene.',
                         'parents': []},
    'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5,
                      'max_positive_instances': None,
                      'max_positive_to_test': None,
                      'max_negative_to_test': None,
                      'max_positive_in_prompt': 50,
                      'max_negative_in_prompt': 20,
                      'max_instances_in_prompt': 100,
                      'test_proportion': 0.1},
    'message': None,
    'attempt': 0,
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': None,
    'num_false_positives': None,
    'num_true_negatives': None,
    'num_false_negatives': None,
    'num_negatives': None,
    'precision': None,
    'recall': None,
    'f1': None,
    'accuracy': None
}


__metadata__ = {   'chemical_class': {   'id': 'CHEBI:35341',
                          'name': 'steroid',
                          'definition': 'Any of naturally occurring compounds '
                                        'and synthetic analogues, based on the '
                                        'cyclopenta[a]phenanthrene carbon '
                                        'skeleton, partially or completely '
                                        'hydrogenated; there are usually '
                                        'methyl groups at C-10 and C-13, and '
                                        'often an alkyl group at C-17. By '
                                        'extension, one or more bond '
                                        'scissions, ring expansions and/or '
                                        'ring contractions of the skeleton may '
                                        'have occurred. Natural steroids are '
                                        'derived biogenetically from squalene '
                                        'which is a triterpene.',
                          'parents': ['CHEBI:18059', 'CHEBI:51958'],
                          'xrefs': ['KEGG:C00377', 'MetaCyc:Steroids'],
                          'all_positive_examples': []},
    'config': None,
    'message': None,
    'sample_true_negatives': [],
    'sample_false_negatives': [   {   'smiles': '[C@@H]1(C[C@H](C(/C(=C/C=C/2\\CCC[C@]3([C@]2(CC[C@]3([H])[C@](CCCC(C)C)([H])C)[H])C)/C1=C)(F)F)O)O',
                                      'name': '4,4-difluoro-1alpha-hydroxyvitamin '
                                              'D3',
                                      'reason': 'Contains 3 rings, less than 4 '
                                                'rings required for steroid '
                                                'backbone'},
                                  {   'smiles': 'C1[C@]2([C@](/C(=C/C=C/3\\C(CC[C@@H](C3)O)=C)/CC1)(CC[C@@]2([C@](C)(CC[C@H](C(C)C)O)O)[H])[H])C',
                                      'name': '(20S,24R)-dihydroxyvitamin D3',
                                      'reason': 'Contains 3 rings, less than 4 '
                                                'rings required for steroid '
                                                'backbone'},
                                  {   'smiles': 'O[C@@H]1C/C(=C\\C=C/2\\[C@]3([C@@]([C@](CC3)([C@@H](CCCC(C)C)C)[H])(CCC2)C)[H])/C(/CC1)=C/C',
                                      'name': '(5E,7E,10E)-(3S)-19-methyl-9,10-seco-5,7,10(19)-cholestatrien-3-ol',
                                      'reason': 'Contains 3 rings, less than 4 '
                                                'rings required for steroid '
                                                'backbone'},
                                  {   'smiles': 'O=C1[C@H](O)[C@]([C@@H]2CCC=3C(C2(C1)C)CC[C@@](C3)(C=C)C)(CO)C',
                                      'name': 'Gifhornenolone B',
                                      'reason': 'Contains 3 rings, less than 4 '
                                                'rings required for steroid '
                                                'backbone'},
                                  {   'smiles': 'O=C(OCCCC)[C@H]1O[C@@H](OC[C@]2([C@@H](O)[C@H](O)C[C@]3([C@H]2CC=C4C3CC[C@@](C4)(C=C)C)C)C)[C@@H](O)[C@@H]([C@@H]1O)O',
                                      'name': 'Virescenoside Z17',
                                      'reason': 'Ring sizes do not match '
                                                'steroid backbone requirements '
                                                '(need at least three '
                                                'six-membered rings and one '
                                                'five-membered ring)'},
                                  {   'smiles': 'O1C(C[C@H]([C@@]2([C@@]3([C@@](CC2)(/C(/CCC3)=C/C=C\\4/C[C@H](O)CCC4=C)[H])C)[H])C)CC(O)(C1O)C',
                                      'name': '(23S,25R)-25-hydroxyvitamin D3 '
                                              '26,23-lactol',
                                      'reason': 'Ring sizes do not match '
                                                'steroid backbone requirements '
                                                '(need at least three '
                                                'six-membered rings and one '
                                                'five-membered ring)'},
                                  {   'smiles': 'O=C1[C@]2(O)[C@@]([C@H]3CCC([C@H]([C@@H]3C1)C=C)=C)(CC[C@H]4[C@@]2(C(=O)O4)C)C',
                                      'name': 'Hawaiinolide A',
                                      'reason': 'Ring sizes do not match '
                                                'steroid backbone requirements '
                                                '(need at least three '
                                                'six-membered rings and one '
                                                'five-membered ring)'},
                                  {   'smiles': 'C1[C@]2([C@](/C(=C/C=C/3\\C(CC[C@@H](C3)O)=C)/CC1)(CC[C@@]2([C@](C)(CC(C(C(C)C)O)O)O)[H])[H])C',
                                      'name': '(20S,23,24)-trihydroxyvitamin '
                                              'D3',
                                      'reason': 'Contains 3 rings, less than 4 '
                                                'rings required for steroid '
                                                'backbone'},
                                  {   'smiles': 'C[C@@H]([C@H]1CC[C@H]2[C@@H]3CC=C4C=CC(=O)[C@]4(C)[C@H]3CC[C@]12C)[C@H]1CC(C)=C(CO)C(=O)O1',
                                      'name': 'withalongolide F',
                                      'reason': 'Fused ring system does not '
                                                'have required ring sizes for '
                                                'steroid backbone'},
                                  {   'smiles': 'O=C1C2=C3N[C@H]4C5=C(C=CC=C5)[C@H]([C@](C3=CC(=C2C(=O)C=6C1=C(O)C=CC6)O)(O)[C@]4(O)[C@H](O)CO)O',
                                      'name': 'Yangpumicin E',
                                      'reason': 'Ring sizes do not match '
                                                'steroid backbone requirements '
                                                '(need at least three '
                                                'six-membered rings and one '
                                                'five-membered ring)'}],
    'attempt': 0,
    'reasoning': 'To classify steroids based on their SMILES strings, we need '
                 'to identify the characteristic steroid backbone, which is a '
                 'tetracyclic fused ring system consisting of three '
                 'six-membered rings and one five-membered ring. Here is the '
                 'Python code using RDKit that implements this classification:',
    'success': True,
    'best': True,
    'error': '',
    'stdout': None,
    'num_true_positives': 973,
    'num_false_positives': 0,
    'num_true_negatives': 0,
    'num_false_negatives': 100,
    'num_negatives': None,
    'precision': 1.0,
    'recall': 0.9068033550792172,
    'f1': 0.9511241446725318,
    'accuracy': 0.9068033550792172,
    'negative_predictive_value': 0.0}