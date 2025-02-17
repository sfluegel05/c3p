"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
#!/usr/bin/env python3
"""
Classifies: Glycosaminoglycan 
(defined as any polysaccharide containing a substantial proportion of aminomonosaccharide residues).

Revised strategy:
1. Check that the molecule is “big” (e.g. molecular weight > 500 Da).
2. Identify candidate sugar rings as non‐aromatic 5‐ or 6‐membered rings that have at least one oxygen atom in the ring
   but no more than 2 total heteroatoms (oxygen and nitrogen) in the ring.
3. For each candidate ring, check whether any ring carbon has a neighbor outside the ring that is a nitrogen (amino substituent).
4. Build connectivity among candidate rings (two rings are “connected” if any atom in one is directly bonded to any atom in the other).
5. Require that the largest connected component of candidate rings has at least three rings and that at least 50% of them are “amino sugar” rings.
6. If these conditions are met, the molecule is classified as a glycosaminoglycan.
 
Note: This is a heuristic algorithm that may still fail in borderline cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is a glycosaminoglycan based on its SMILES string.
    
    The approach is to (1) ensure the molecule is big enough (MW > 500),
    (2) identify candidate sugar rings using a less strict rule (non‐aromatic, 5‐ or 6‐membered, with 1–2 ring heteroatoms),
    (3) build a connectivity graph of such rings (requiring at least 3 to call it a polysaccharide),
    and (4) check that at least 50% of those rings have an external amino substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as a glycosaminoglycan, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Step 0: Check if molecular weight is high enough.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a polysaccharide."

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_rings = []      # each candidate is a set of atom indices
    ring_amino_flags = []     # flag for each candidate: True if ring shows external amino substitution

    # Step 1: Identify candidate sugar rings.
    # Accept rings of size 5 or 6 that are non‐aromatic and have at least one oxygen in the ring.
    # Also allow at most 2 heteroatoms (oxygen and nitrogen) in the ring.
    for ring in atom_rings:
        if len(ring) not in (5, 6):
            continue
        # Get atoms corresponding to indices in the ring.
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Skip if any atom in the ring is aromatic.
        if any(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue
        # Count heteroatoms in the ring (only consider O and N).
        hetero_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() in (7, 8))
        # Require at least one oxygen.
        oxygen_count = sum(1 for atom in atoms_in_ring if atom.GetAtomicNum() == 8)
        if oxygen_count < 1:
            continue
        # Accept rings that have no more than 2 heteroatoms (to avoid rings with extra heteroatoms not typical for sugars).
        if hetero_count > 2:
            continue
        
        candidate_rings.append(set(ring))
        
        # Step 2: Check for external amino substituents.
        # For each carbon in the ring, see if it has a neighbor (outside the ring) that is nitrogen.
        amino_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Typically we look at carbons. (Many sugars are mainly carbons.)
            if atom.GetAtomicNum() != 6:
                continue
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue  # skip atoms that are part of the ring
                if neighbor.GetAtomicNum() == 7:
                    amino_found = True
                    break
            if amino_found:
                break
        ring_amino_flags.append(amino_found)

    total_candidates = len(candidate_rings)
    if total_candidates == 0:
        return False, "No candidate sugar rings found."
    
    # Step 3: Build connectivity among candidate rings.
    # Two candidate rings are considered connected if any atom in one is bonded to any atom in the other.
    n = total_candidates
    connections = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            # Check each pair of atoms (one in ring i, one in ring j) to see if there is a bond.
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
    
    # Find connected components among candidate rings using depth-first search.
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
    
    # Get the largest connected component.
    largest_comp = max(components, key=len) if components else set()
    comp_size = len(largest_comp)
    if comp_size < 3:
        return False, f"Not a polysaccharide: largest connected candidate sugar network has only {comp_size} ring(s)."
    
    # Step 4: Within the largest connected network, count how many rings have an external amino substituent.
    comp_amino = sum(1 for i in largest_comp if ring_amino_flags[i])
    ratio = comp_amino / comp_size

    if ratio >= 0.5:
        return True, (f"Polysaccharide with {comp_size} connected sugar ring(s), "
                      f"of which {comp_amino} ({ratio*100:.1f}%) show external amino substitution.")
    else:
        return False, (f"Candidate polysaccharide with {comp_size} connected sugar ring(s) "
                       f"but only {comp_amino} ({ratio*100:.1f}%) appear to be amino sugars.")

# Example (for testing purposes):
if __name__ == "__main__":
    # Example SMILES (one of the given glycosaminoglycan examples)
    test_smiles = "O=C1O[C@@H](C=C[C@@H](CC(C=O)C=2C3=C([C@@]4([C@H](C=C1C)C(=O)NC4=CC3=O)O)C=C(C)C2O)CC)[C@@H](O)C=C(C)C"  # Divergolide D
    result, reason = is_glycosaminoglycan(test_smiles)
    print("Result:", result)
    print("Reason:", reason)