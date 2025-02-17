"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
#!/usr/bin/env python3
"""
Classifies: Glycosaminoglycan
Definition: Any polysaccharide containing a substantial proportion of aminomonosaccharide residues.

Improved strategy:
1. Check that the molecule is “big enough” (here, MW > 300 Da).
2. Identify candidate sugar rings as non‐aromatic, 5‐ or 6‐membered rings that:
   a. Have at least one ring oxygen and at most 2 heteroatoms (O and N).
   b. Contain only single bonds in the ring.
   c. Display at least 1 external substituent of hydroxyl (O with H) or amino (N with H) type.
3. For each candidate ring, flag it if at least one ring carbon has an external amino substituent.
4. Build connectivity among candidate rings (two rings are “connected” if any atom of one is bonded to any atom of the other).
5. Require that the largest connected component of candidate rings has at least three rings and that at least 30% of them have the amino flag.
6. If these conditions are met, classify the molecule as a glycosaminoglycan.

Note: This heuristic algorithm is imperfect.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines if a molecule is classified as a glycosaminoglycan
    (i.e. a polysaccharide containing a substantial proportion
    of aminomonosaccharide residues) based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as glycosaminoglycan, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1: Basic check on molecular weight (using a lowered threshold)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a polysaccharide."
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_rings = []   # List of sets containing atom indices for candidate rings.
    ring_amino_flags = []  # Parallel list to flag if candidate ring has an amino substitution.

    # Iterate over all rings
    for ring in atom_rings:
        # Consider only rings of size 5 or 6
        if len(ring) not in (5, 6):
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Exclude rings containing aromatic atoms
        if any(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue
        
        # Check that all bonds in the ring are single bonds
        all_single = True
        for i in range(len(ring)):
            a_idx = ring[i]
            b_idx = ring[(i+1) % len(ring)]
            bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                all_single = False
                break
        if not all_single:
            continue
        
        # Count heteroatoms (only oxygen and nitrogen)
        hetero_atoms = [atom for atom in atoms_in_ring if atom.GetAtomicNum() in (7,8)]
        oxygen_atoms = [atom for atom in atoms_in_ring if atom.GetAtomicNum() == 8]
        if len(oxygen_atoms) < 1:
            continue
        # Allow at most 2 heteroatoms in the ring (to allow one extra nitrogen in aminomonosaccharides)
        if len(hetero_atoms) > 2:
            continue

        # Count external substituents (neighbors outside the ring that are O or N with at least one H)
        substituent_count = 0
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                if neighbor.GetAtomicNum() in (7,8) and neighbor.GetTotalNumHs() > 0:
                    substituent_count += 1
        if substituent_count < 1:
            continue

        candidate_rings.append(set(ring))
        
        # Flag the candidate ring if any ring carbon has an external amino substituent.
        amino_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # normally sugar rings contain C. If the atom is carbon, check its neighbors.
            if atom.GetAtomicNum() != 6:
                continue
            for neighbor in atom.GetNeighbors():
                if neighbor.GetIdx() in ring:
                    continue
                # A candidate amino substituent is a nitrogen with at least one attached hydrogen.
                if neighbor.GetAtomicNum() == 7 and neighbor.GetTotalNumHs() > 0:
                    amino_found = True
                    break
            if amino_found:
                break
        ring_amino_flags.append(amino_found)
    
    total_candidates = len(candidate_rings)
    if total_candidates == 0:
        return False, "No candidate sugar rings found."
    
    # Build connectivity between candidate rings.
    # Two rings are connected if any atom from one ring is bonded to any atom of the other.
    n = total_candidates
    connections = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
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
    
    # Identify connected components using DFS.
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
    
    # Count how many rings in the connected component show an external amino substituent.
    comp_amino = sum(1 for i in largest_comp if ring_amino_flags[i])
    ratio = comp_amino / comp_size

    # Lower threshold to 30% for classification.
    if ratio >= 0.3:
        return True, (f"Polysaccharide with {comp_size} connected sugar ring(s), "
                      f"of which {comp_amino} ({ratio*100:.1f}%) show external amino substitution.")
    else:
        return False, (f"Candidate polysaccharide with {comp_size} connected sugar ring(s) "
                       f"but only {comp_amino} ({ratio*100:.1f}%) appear to be amino sugars.")

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Try one of the provided SMILES strings.
    test_smiles = "O=C1O[C@@H](C=C[C@@H](CC(C=O)C=2C3=C([C@@]4([C@H](C=C1C)C(=O)NC4=CC3=O)O)C=C(C)C2O)CC)[C@@H](O)C=C(C)C"  # Example: Divergolide D
    result, reason = is_glycosaminoglycan(test_smiles)
    print("Result:", result)
    print("Reason:", reason)