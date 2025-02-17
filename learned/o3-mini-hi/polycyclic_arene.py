"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: Polycyclic Arene – defined as a polycyclic aromatic hydrocarbon (PAH).

This improved version revises the previous algorithm by:
  • Lowering the overall heavy atom threshold to 10 (so that small PAHs like naphthalene are included).
  • Computing coverage over the molecule’s aromatic carbon atoms rather than all heavy atoms.
  • Lowering the minimum coverage threshold for the fused aromatic candidate system (from 0.55 to 0.50).
  • Requiring that a fused system is made of at least 2 candidate rings and at least 7 atoms.
  
A molecule is classified as a PAH if there is a connected (i.e. fused) group of rings, where each candidate ring is composed solely of aromatic carbons,
and the union of these rings covers at least 50% of the molecule’s aromatic carbon atoms.
"""

from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon, PAH)
    based on its SMILES string.
    
    The algorithm:
      - Parses the SMILES string.
      - Computes the set of aromatic carbon atoms in the molecule.
      - Obtains ring information from the molecule.
      - Flags each ring as a candidate if every atom in the ring is aromatic and is carbon.
      - Constructs a connectivity graph among rings, where two rings are fused if they share at least 2 atoms.
      - For each connected component of candidate rings, if there are at least 2 candidate rings and
        their union (the fused candidate system) has at least 7 atoms, then we compute its “coverage” – the fraction
        of aromatic carbon atoms of the molecule that lie in those fused rings. If this coverage is at least 0.50,
        the molecule is classified as a PAH.
        
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        (bool, str): Tuple. First element is True if the molecule is classified as a PAH, False otherwise;
                     second element gives the rationale.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get aromatic carbons (only consider atoms that are aromatic and carbon)
    aromatic_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic()]
    total_aromatic_carb = len(aromatic_carbons)
    if total_aromatic_carb < 4:  # too few aromatic carbons; most PAHs will be larger than naphthalene (10 aromatic C's)
        return False, f"Total aromatic carbon count ({total_aromatic_carb}) is too low for a PAH"
    
    # Also check overall heavy atom count; lower threshold to include naphthalene
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    total_heavy = len(heavy_atoms)
    if total_heavy < 10:
        return False, f"Total heavy atom count ({total_heavy}) is below the minimum required for a PAH"
    
    # Get the ring information from the molecule.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the molecule"
    
    # Candidate ring: each ring must be made entirely of aromatic carbons.
    def is_candidate_ring(ring):
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not (atom.GetIsAromatic() and atom.GetAtomicNum() == 6):
                return False
        return True

    candidate_flags = [is_candidate_ring(ring) for ring in rings]
    ring_atom_sets = [set(ring) for ring in rings]
    
    # Build connectivity graph among rings
    # Two rings are considered fused if they share at least 2 atoms.
    n = len(rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if len(ring_atom_sets[i].intersection(ring_atom_sets[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Tunable thresholds:
    min_component_candidate_rings = 2    # at least 2 candidate rings in a fused aromatic system
    min_candidate_atoms = 7              # fused system must have at least 7 atoms
    coverage_threshold = 0.50            # fused system must cover at least 50% of aromatic carbons in the molecule
    
    # Search for connected components using a depth-first search
    visited = set()
    for i in range(n):
        if i in visited:
            continue
        # Get the connected component via DFS
        stack = [i]
        component = set()
        while stack:
            current = stack.pop()
            if current not in component:
                component.add(current)
                stack.extend(ring_graph[current] - component)
        visited |= component
        
        # Count candidate rings in this component.
        candidate_count = sum(1 for idx in component if candidate_flags[idx])
        if candidate_count < min_component_candidate_rings:
            continue
        
        # Build union of atoms only for candidate rings in this component.
        candidate_atoms = set()
        for idx in component:
            if candidate_flags[idx]:
                candidate_atoms |= ring_atom_sets[idx]
        if len(candidate_atoms) < min_candidate_atoms:
            continue
        
        # Now compute coverage as ratio of candidate aromatic carbons (from fused candidate rings)
        # to the molecule’s total aromatic carbons.
        coverage = len(candidate_atoms.intersection({atom.GetIdx() for atom in aromatic_carbons})) / total_aromatic_carb
        
        if coverage < coverage_threshold:
            # Although a fused candidate system was detected, it covers a low fraction of aromatic carbons.
            return False, (f"Found a fused aromatic candidate system with {candidate_count} candidate rings, "
                           f"but it only covers {coverage:.2f} of the molecule's aromatic carbons")
        
        # If criteria are met, classify as PAH.
        return True, (f"Found a fused aromatic hydrocarbon system with {candidate_count} candidate rings "
                      f"covering {coverage:.2f} of the aromatic carbons (candidate system atoms: {len(candidate_atoms)})")
        
    return False, "No fused aromatic hydrocarbon system detected"


# Example testing code (can be removed when integrating)
if __name__ == '__main__':
    test_smiles = [
        "c1ccc2c(c1)ccc1cc3c(ccc4ccccc34)cc21",  # dibenz[a,h]anthracene (expected True)
        "c1ccc2cc3ccccc3cc2c1",                  # anthracene (expected True)
        "Oc1cc2ccc3cc4ccccc4c4ccc(c1)c2c34",      # 2-Hydroxybenzo[a]pyrene (expected True)
        "Cn1c2ccccc2c(=O)c2c(O)cc3occc3c12",       # Furofoline I (expected False)
        "Clc1cc(Cl)c2[C@@H](c3ccccc3CCc2c1)n1ccnc1",  # (R)-eberconazole (expected False)
        "c1ccc2c(c1)cc1C=Cc3cccc2c13",            # acephenanthrylene (expected True)
        "c1ccc2ccccc2c1",                        # naphthalene (expected True)        
        "Cc1cc(O)c2c3c1c1c(CO)cc(O)c4c1c1c5c(c(O)cc(O)c5c4=O)c4c(O)cc(O)c(c4c31)c2=O"  # pseudohypericin (expected True)
    ]
    
    for smi in test_smiles:
        result, reason = is_polycyclic_arene(smi)
        print(smi, "=>", result, ":", reason)