"""
Classifies: CHEBI:18085 glycosaminoglycan
"""
#!/usr/bin/env python3
"""
Classifies: Glycosaminoglycan
Definition: Any polysaccharide containing a substantial proportion of aminomonosaccharide residues.

This updated algorithm works as follows:
1. Check that the molecule is “big enough” (here, MW > 300 Da).
2. Identify candidate sugar rings as non‐aromatic 5‐ or 6‐membered rings that:
   a. Contain at least one ring oxygen.
   b. Contain at most three heteroatoms (O or N) in the ring.
   c. Have at least 80% of the atoms in the ring with sp3 hybridization.
   d. Possess at least one external substituent that is hydroxyl (O–H) or amino (N–H).
3. For each candidate ring, flag it if any ring carbon carries an external amino substituent.
4. Build connectivity among candidate rings (two rings are “connected” if any atom in one is bonded to any atom in the other).
5. Require that the largest connected component of candidate rings has at least three rings and that at least 30% of them show the amino flag.
6. If so, classify as glycosaminoglycan.

Note: This heuristic algorithm is imperfect.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycosaminoglycan(smiles: str):
    """
    Determines whether a molecule qualifies as a glycosaminoglycan
    (i.e. a polysaccharide with a substantial proportion of aminomonosaccharide residues)
    based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if classified as glycosaminoglycan, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Step 1: Basic molecular weight check
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) to be a polysaccharide."
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    
    candidate_rings = []   # Each candidate ring as a set of atom indices.
    ring_amino_flags = []  # Parallel list that flags if candidate ring has a carbon externally substituted with an amino group.
    
    # Iterate over all rings
    for ring in atom_rings:
        # Consider only 5- or 6-membered rings.
        if len(ring) not in (5, 6):
            continue
        atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # Exclude aromatic rings.
        if any(atom.GetIsAromatic() for atom in atoms_in_ring):
            continue
        
        # Criterion A: The ring should contain at least one oxygen.
        ring_oxygens = [atom for atom in atoms_in_ring if atom.GetAtomicNum() == 8]
        if len(ring_oxygens) < 1:
            continue
        
        # Criterion B: Allow at most three heteroatoms (O or N) in the ring.
        hetero_atoms = [atom for atom in atoms_in_ring if atom.GetAtomicNum() in (7,8)]
        if len(hetero_atoms) > 3:
            continue
        
        # Criterion C: Require that at least 80% of atoms in the ring are sp3-hybridized.
        sp3_count = sum(1 for atom in atoms_in_ring if atom.GetHybridization() == Chem.HybridizationType.SP3)
        if sp3_count / len(ring) < 0.8:
            continue
        
        # Criterion D: Count external substituents that are either hydroxyl (O with at least one H)
        # or amino (N with at least one H) groups.
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
        
        # Passed candidate ring tests: add candidate.
        candidate_rings.append(set(ring))
        
        # Flag: For candidate ring, if any ring carbon (atomic number 6) has an external amino (-NH2)
        # substituent, flag this ring.
        amino_found = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
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
    
    # Build connectivity network among candidate rings.
    # Two rings are connected if any atom from one is bonded to any atom from the other.
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
    
    # Count rings in the component that show an amino substitution.
    comp_amino = sum(1 for i in largest_comp if ring_amino_flags[i])
    ratio = comp_amino / comp_size
    
    # Classification: at least 30% of the connected rings should be amino sugars.
    if ratio >= 0.3:
        return True, (f"Polysaccharide with {comp_size} connected sugar ring(s), "
                      f"of which {comp_amino} ({ratio*100:.1f}%) show external amino substitution.")
    else:
        return False, (f"Candidate polysaccharide with {comp_size} connected sugar ring(s) "
                       f"but only {comp_amino} ({ratio*100:.1f}%) appear to be amino sugars.")

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Test with one provided example (Divergolide D)
    test_smiles = "O=C1O[C@@H](C=C[C@@H](CC(C=O)C=2C3=C([C@@]4([C@H](C=C1C)C(=O)NC4=CC3=O)O)C=C(C)C2O)CC)[C@@H](O)C=C(C)C"
    result, reason = is_glycosaminoglycan(test_smiles)
    print("Result:", result)
    print("Reason:", reason)