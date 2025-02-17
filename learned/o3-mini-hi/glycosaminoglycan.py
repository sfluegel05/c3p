"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
#!/usr/bin/env python3
"""
Classifies: Glycosaminoglycan (defined as any polysaccharide containing 
a substantial proportion of aminomonosaccharide residues)

This improved version uses a multi‐step heuristic:
  1. Identify candidate sugar rings as 5‐ or 6‐membered rings that have exactly one oxygen atom in the ring.
  2. Build a connectivity graph among candidate rings based on shared bonds outside the ring.
  3. In order to be considered a polysaccharide, the largest connected set of candidate rings must have at least 3 rings.
  4. For each candidate ring in that connected network, check whether any ring carbon has a substituent nitrogen (not in the ring).
  5. Declare the molecule a glycosaminoglycan if at least 50% of the connected sugar rings are “amino sugars.”
  
Note: This is a heuristic classifier and may still fail in borderline cases.
"""

from rdkit import Chem

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    The approach is to (1) identify candidate sugar rings and (2) verify 
    that they form a connected polysaccharide network and (3) that a large portion 
    of those rings have an external nitrogen substituent.

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a glycosaminoglycan, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of atom-index tuples

    candidate_rings = []       # list of rings (each is set of atom indices)
    ring_amino_flags = []      # corresponding flag: True if ring has an external amino substituent

    # Step 1: Identify candidate sugar rings.
    # We consider rings of size 5 or 6 having exactly one oxygen atom in the ring.
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue
        # Count the number of oxygen atoms in the ring.
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        oxygen_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 8)
        if oxygen_count != 1:
            continue
        candidate_rings.append(set(ring))
        # Check for amino substitution: for each carbon in the ring, 
        # see if there is any neighbor outside the ring that is a nitrogen.
        amino_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # only check carbons
                continue
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in ring:
                    continue   # skip atoms within the ring itself
                if nb.GetAtomicNum() == 7:
                    amino_found = True
                    break
            if amino_found:
                break
        ring_amino_flags.append(amino_found)

    total_candidates = len(candidate_rings)
    if total_candidates == 0:
        return False, "No candidate sugar rings found."

    # Step 2: Build connectivity among candidate rings.
    # Two candidate rings are considered connected if they share any bond through substituents.
    # Here, if any atom in ring A is directly bonded (in the molecule) to an atom in ring B 
    # (and those atoms are not the same, i.e. each belongs to its candidate ring), add an edge.
    n = total_candidates
    connections = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            # Check if ring i and ring j are connected by any bond.
            # We say they are connected if any atom in candidate_rings[i] is bonded
            # to any atom in candidate_rings[j].
            if any(mol.GetBondBetweenAtoms(a, b) is not None for a in candidate_rings[i] for b in candidate_rings[j]):
                connections[i].add(j)
                connections[j].add(i)

    # Find connected components using a depth-first search.
    visited = set()
    components = []
    for i in range(n):
        if i in visited:
            continue
        comp = set()
        stack = [i]
        while stack:
            curr = stack.pop()
            if curr in comp:
                continue
            comp.add(curr)
            visited.add(curr)
            stack.extend(connections[curr] - comp)
        components.append(comp)

    # Get the largest connected component
    largest_comp = max(components, key=len) if components else set()
    comp_size = len(largest_comp)
    if comp_size < 3:
        return False, f"Not a polysaccharide: largest connected candidate sugar network has only {comp_size} ring(s)."
    
    # Within the largest connected network, count the number of candidate rings and how many are amino sugars.
    comp_total = comp_size
    comp_amino = sum(1 for i in largest_comp if ring_amino_flags[i])
    ratio = comp_amino / comp_total

    if ratio >= 0.5:
        return True, (f"Polysaccharide with {comp_total} connected sugar ring(s), "
                      f"of which {comp_amino} ({ratio*100:.1f}%) contain an amino substituent.")
    else:
        return False, (f"Candidate polysaccharide with {comp_total} connected sugar ring(s) "
                       f"but only {comp_amino} ({ratio*100:.1f}%) are amino sugars.")

# Example usage:
if __name__ == "__main__":
    # For example, test one of the provided SMILES strings.
    test_smiles = "O=C1O[C@@H](C=C[C@@H](CC[C@]23O[C@H](C=4C(=C(NC(CC=C1C)=O)C=C(O)C4)O2)[C@@H](C)C(C3)=O)CC)[C@@H](O)C=C(C)C"
    result, reason = is_glycosaminoglycan(test_smiles)
    print("Result:", result)
    print("Reason:", reason)