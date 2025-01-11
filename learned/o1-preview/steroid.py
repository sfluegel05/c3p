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