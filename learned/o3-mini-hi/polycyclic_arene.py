"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: Polycyclic Arene, defined as a polycyclic aromatic hydrocarbon (PAH).

In this version the algorithm improves on the previous attempt by:
  1. Lowering the coverage threshold of candidate ring atoms (from 0.60 to 0.55).
  2. Allowing a smaller fused candidate system (at least 7 atoms).
  3. Requiring a minimum total heavy atom count (≥12) so that very small systems are rejected.
  4. Reporting the number of fused candidate rings.
  
A molecule is classified as a PAH if there is a connected (i.e. fused)
system of candidate rings (where every atom in each ring is aromatic carbon)
meeting the above criteria.
"""

from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon, PAH)
    based on its SMILES string.

    The algorithm performs the following steps:
      - Parse the SMILES string and extract ring information.
      - Define candidate rings as those rings in which each atom is aromatic and is carbon.
      - Build a connectivity graph linking rings that share at least 2 atoms.
      - For each connected component of candidate rings, require:
            • At least 2 candidate rings,
            • Their union (fused candidate system) contains at least 7 atoms,
            • The fused system covers at least 55% of the molecule's heavy atoms,
            • And the overall molecule has at least 12 heavy atoms.
      - If any connected fused candidate system satisfies these, classify as PAH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        (bool, str): Tuple. First element is True if the molecule is classified as a PAH,
                     False otherwise; second element gives the rationale.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count heavy (non-H) atoms.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    total_heavy = len(heavy_atoms)
    if total_heavy < 12:
        return False, f"Total heavy atom count ({total_heavy}) is below the minimum required for a PAH"
    
    # Obtain ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the molecule"
    
    # Helper function: determine if a given ring is a candidate aromatic ring.
    # It must have every atom aromatic and be a carbon.
    def is_candidate_ring(ring):
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic() or atom.GetAtomicNum() != 6:
                return False
        return True

    # Flag each ring whether it qualifies as a candidate ring.
    candidate_flags = [is_candidate_ring(ring) for ring in rings]

    # Build a connectivity graph among all rings.
    # Two rings are considered fused if they share at least 2 atoms.
    n = len(rings)
    ring_graph = {i: set() for i in range(n)}
    ring_atom_sets = [set(r) for r in rings]
    for i in range(n):
        for j in range(i + 1, n):
            if len(ring_atom_sets[i].intersection(ring_atom_sets[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)

    # We now search for connected components from the rings.
    visited = set()
    # Tunable thresholds:
    coverage_threshold = 0.55  # fraction of molecule heavy atoms covered by fused PAH core
    min_candidate_atoms = 7    # minimum number of atoms in the fused candidate ring system
    
    for i in range(n):
        if i in visited:
            continue
        # Depth-first traversal to get a connected component (fused system)
        stack = [i]
        component = set()
        while stack:
            current = stack.pop()
            if current not in component:
                component.add(current)
                stack.extend(ring_graph[current] - component)
        visited |= component

        # Count how many rings in this component are candidate rings.
        candidate_count = sum(1 for idx in component if candidate_flags[idx])
        if candidate_count < 2:
            continue  # Need at least 2 candidate rings to be considered a fused aromatic system.
        
        # Build the union of atoms from the candidate rings in the component.
        candidate_atoms = set()
        for idx in component:
            if candidate_flags[idx]:
                candidate_atoms |= ring_atom_sets[idx]
        if len(candidate_atoms) < min_candidate_atoms:
            continue  # fused candidate system is too small
        
        coverage = len(candidate_atoms) / total_heavy
        if coverage < coverage_threshold:
            # Although we found a fused candidate system, it only covers a small fraction
            # of the heavy atoms in the molecule.
            return False, (f"Found a fused aromatic candidate system with {candidate_count} candidate rings, "
                           f"but it only covers {coverage:.2f} of the molecule's heavy atoms")
        
        # If all criteria are met then return success.
        return True, (f"Found a fused aromatic hydrocarbon system with {candidate_count} candidate rings "
                      f"covering {coverage:.2f} of the heavy atoms (candidate system atoms: {len(candidate_atoms)})")
    
    return False, "No fused aromatic hydrocarbon system detected"


# Example testing code (can be removed when integrating)

if __name__ == '__main__':
    test_smiles = [
        "c1ccc2c(c1)ccc1cc3c(ccc4ccccc34)cc21",  # dibenz[a,h]anthracene (expected True)
        "c1ccc2cc3ccccc3cc2c1",                  # anthracene (expected True)
        "Oc1cc2ccc3cc4ccccc4c4ccc(c1)c2c34",      # 2-Hydroxybenzo[a]pyrene (expected True)
        "N(C)C1C=2C(N(C3=C(C1)C=CC=C3)C)=CC=CC2",  # metapramine (expected False)
        "Cn1c2ccccc2c(=O)c2c(O)cc3occc3c12",       # Furofoline I (expected False)
        "Cc1cc(O)c2c3c1c1c(CO)cc(O)c4c1c1c5c(c(O)cc(O)c5c4=O)c4c(O)cc(O)c(c4c31)c2=O"  # pseudohypericin (expected True)
    ]
    
    for smi in test_smiles:
        result, reason = is_polycyclic_arene(smi)
        print(smi, "=>", result, ":", reason)