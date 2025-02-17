"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
#!/usr/bin/env python3
"""
Classifies: Glycosaminoglycan 
(defined as any polysaccharide containing a substantial proportion of aminomonosaccharide residues).

Improved strategy:
1. Check that the molecule is “big” (e.g. molecular weight > 500 Da).
2. Identify candidate sugar rings as non‐aromatic, 5‐ or 6‐membered rings that:
   a. Have at least one ring oxygen and at most 2 total heteroatoms (O and N).
   b. Contain only single bonds (to capture the typical saturated sugar rings).
   c. Show at least 2 external substituents that are either hydroxyl (O with at least one hydrogen) or amino groups.
3. For each candidate ring, flag it if at least one ring carbon has an external amino substituent 
   (a neighboring N atom that carries at least one hydrogen).
4. Build connectivity among candidate rings (two rings are “connected” if any atom in one is bonded to any atom in the other).
5. Require that the largest connected component of candidate rings has at least three rings and that at least 50% of them are “amino sugar” rings.
6. If these conditions are met, classify the molecule as a glycosaminoglycan.
Note: This heuristic algorithm is imperfect.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    
    The approach is to (1) ensure the molecule is big enough (MW > 500 Da),
    (2) identify candidate sugar rings that are non‐aromatic 5‐ or 6‐membered rings with only single bonds,
        have at least one oxygen, no more than 2 heteroatoms (N or O),
        and that display at least 2 external substituents that are hydroxyl or amino groups,
    (3) flag rings that show an external amino substituent (N with at least one H),
    (4) build a connectivity graph of such rings (requiring at least 3 to call it a polysaccharide),
    and (5) check that at least 50% of those rings have the amino flag.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a glycosaminoglycan, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Step 0: Basic check on molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a polysaccharide."

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_rings = []      # Each candidate is a set of atom indices
    ring_amino_flags = []     # For each candidate, True if ring shows external amino substitution
    
    # Iterate over all rings
    for ring in atom_rings:
        # Only consider rings of size 5 or 6
        if len(ring) not in (5, 6):
            continue
        # Get atoms in the ring
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Exclude rings with any aromatic atom
        if any(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue
        
        # Check that all bonds in the ring are single (sugar rings are generally saturated)
        bonds_in_ring = []
        for i in range(len(ring)):
            a_idx = ring[i]
            b_idx = ring[(i+1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                bonds_in_ring = None
                break
            bonds_in_ring.append(bond)
        if bonds_in_ring is None:
            continue  # Skip rings with non-single bonds

        # Count heteroatoms (only consider oxygen and nitrogen)
        hetero_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() in (7, 8))
        # At least one oxygen must be present
        oxygen_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 8)
        if oxygen_count < 1:
            continue
        if hetero_count > 2:
            continue

        # Count external substituents on ring atoms that look like hydroxyl or amino groups.
        substituent_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                # Consider only neighbors outside the ring.
                if neighbor.GetIdx() in ring:
                    continue
                # We'll assume a substituent is hydroxyl or amino if it is O or N and has at least one hydrogen.
                if neighbor.GetAtomicNum() in (7, 8) and neighbor.GetTotalNumHs() > 0:
                    substituent_count += 1
        # Require at least 2 such substituents to favor sugar-like topology.
        if substituent_count < 2:
            continue

        # Add the ring as a candidate.
        candidate_rings.append(set(ring))
        
        # Check for external amino substituents:
        amino_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look typically at carbon atoms in the ring.
            if atom.GetAtomicNum() != 6:
                continue
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                if neighbor.GetAtomicNum() == 7 and neighbor.GetTotalNumHs() > 0:
                    amino_found = True
                    break
            if amino_found:
                break
        ring_amino_flags.append(amino_found)
    
    total_candidates = len(candidate_rings)
    if total_candidates == 0:
        return False, "No candidate sugar rings found."
    
    # Build connectivity among candidate rings.
    # Two rings are considered connected if any atom in one is bonded to any atom in the other.
    n = total_candidates
    connections = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            connected = False
            for a in candidate_rings[i]:
                for b in candidate_rings[j]:
                    if mol.GetBondBetweenAtoms(a, b) is not None:
                        connected = True
                        break
                if connected:
                    break
            if connected:
                connections[i].add(j)
                connections[j].add(i)
    
    # Find connected components among candidate rings using DFS.
    visited = set()
    components = []
    for i in range(n):
        if i in visited:
            continue
        comp = set()
        stack = [i]
        while stack:
            cur = stack.pop()
            if cur in comp:
                continue
            comp.add(cur)
            visited.add(cur)
            stack.extend(connections[cur] - comp)
        components.append(comp)
    
    # Select the largest connected component.
    largest_comp = max(components, key=len) if components else set()
    comp_size = len(largest_comp)
    if comp_size < 3:
        return False, f"Not a polysaccharide: largest connected candidate sugar network has only {comp_size} ring(s)."
    
    # Within the largest connected network, count rings with an external amino substituent.
    comp_amino = sum(1 for i in largest_comp if ring_amino_flags[i])
    ratio = comp_amino / comp_size

    if ratio >= 0.5:
        return True, (f"Polysaccharide with {comp_size} connected sugar ring(s), "
                      f"of which {comp_amino} ({ratio*100:.1f}%) show external amino substitution.")
    else:
        return False, (f"Candidate polysaccharide with {comp_size} connected sugar ring(s) "
                       f"but only {comp_amino} ({ratio*100:.1f}%) appear to be amino sugars.")

# Example usage (for testing purposes):
if __name__ == "__main__":
    # You can test with one of the provided SMILES strings.
    test_smiles = "O=C1O[C@@H](C=C[C@@H](CC(C=O)C=2C3=C([C@@]4([C@H](C=C1C)C(=O)NC4=CC3=O)O)C=C(C)C2O)CC)[C@@H](O)C=C(C)C"  # Example: Divergolide D
    result, reason = is_glycosaminoglycan(test_smiles)
    print("Result:", result)
    print("Reason:", reason)