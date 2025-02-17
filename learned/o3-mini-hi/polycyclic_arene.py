"""
Classifies: CHEBI:33848 polycyclic arene
"""
"""
Classifies: Polycyclic Arene, defined as a polycyclic aromatic hydrocarbon (PAH).

An entity is classified as a PAH if it contains a fused system of aromatic rings that are
mostly made up of carbons. In our algorithm we:
  1. Parse the SMILES and extract all rings.
  2. Flag candidate rings as those in which every atom is aromatic and is a carbon.
  3. Build a “ring connectivity graph” by linking rings that share at least 2 atoms.
  4. Identify connected components of candidate rings.
  5. For each such component, require:
       • at least 2 candidate rings,
       • a fused candidate system (union of candidate ring atoms) of size ≥10,
       • the candidate system covers a significant fraction (≥ 0.60) of the molecule’s heavy atoms.
  6. If these criteria are met we classify the molecule as a PAH.
  
This method is intended to avoid false positives (small aromatic substituents) while still
capturing PAHs that have heteroatom substituents on an otherwise fused aromatic core.
"""

from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (polycyclic aromatic hydrocarbon, PAH)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if the molecule is classified as a PAH,
                     and the second gives a reason.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the molecule"
    
    # Helper: Check if a ring is a candidate aromatic (core) ring.
    # In our definition a candidate ring should have every atom aromatic and be carbon.
    def is_candidate_ring(ring):
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic() or (atom.GetAtomicNum() != 6):
                return False
        return True

    # Mark each ring as a candidate (True/False).
    candidate_flags = [is_candidate_ring(ring) for ring in rings]
    
    # Build a connectivity graph among rings.
    # Two rings are considered fused if they share at least 2 atoms.
    n = len(rings)
    ring_graph = {i: set() for i in range(n)}
    ring_atom_sets = [set(r) for r in rings]
    for i in range(n):
        for j in range(i+1, n):
            if len(ring_atom_sets[i].intersection(ring_atom_sets[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components in the ring graph.
    visited = set()
    for i in range(n):
        if i in visited:
            continue
        # Do a depth-first search for the connected component.
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
        if candidate_count >= 2:
            # Build the union of atoms in all candidate rings in the component.
            candidate_atoms = set()
            for idx in component:
                if candidate_flags[idx]:
                    candidate_atoms |= ring_atom_sets[idx]
            if len(candidate_atoms) < 10:
                # Too small a candidate system to be a true PAH.
                continue

            # Now check that the fused candidate system covers a large fraction of the molecule.
            # Calculate the total number of heavy atoms (non-hydrogen).
            heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
            total_heavy = len(heavy_atoms)
            if total_heavy == 0:
                return False, "No heavy atoms found in the molecule"
            coverage = len(candidate_atoms) / total_heavy
            
            # We require the candidate fused system to account for at least 60% of heavy atoms.
            if coverage < 0.60:
                return False, (f"Found a fused aromatic candidate system with {candidate_count} candidate rings, "
                               f"but it only covers {coverage:.2f} of the molecule's heavy atoms")
            
            return True, f"Found a fused aromatic hydrocarbon system with {candidate_count} candidate rings"

    return False, "No fused aromatic hydrocarbon system detected"


# Example testing code (these lines may be removed or commented out when integrating)
if __name__ == '__main__':
    test_smiles = [
        "c1ccc2c(c1)ccc1cc3c(ccc4ccccc34)cc21",  # dibenz[a,h]anthracene (expected True)
        "c1ccc2cc3ccccc3cc2c1",                  # anthracene (expected True)
        "Oc1cc2ccc3cc4ccccc4c4ccc(c1)c2c34",      # 2-Hydroxybenzo[a]pyrene (expected True)
        "N(C)C1C=2C(N(C3=C(C1)C=CC=C3)C)=CC=CC2",  # metapramine (expected False)
        "Cc1cc(O)c2c3c1c1c(CO)cc(O)c4c1c1c5c(c(O)cc(O)c5c4=O)c4c(O)cc(O)c(c4c31)c2=O"  # pseudohypericin (expected True)
    ]
    
    for smi in test_smiles:
        result, reason = is_polycyclic_arene(smi)
        print(smi, "=>", result, ":", reason)