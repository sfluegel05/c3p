"""
Classifies: CHEBI:26125 phytosterols
"""
#!/usr/bin/env python
"""
Classifies: Phytosterols – Sterols similar to cholesterol that occur in plants.
They are defined here as molecules that contain a fused steroid nucleus 
(a tetracyclic system of 4 rings – one five‐membered and three six‐membered rings,
all made of carbon and containing roughly 17 carbons in the nucleus),
a hydroxyl group (–OH) typically attached to that core,
and an aliphatic side chain attached to the nucleus that is at least 4 carbons long.
The overall molecular weight should be in the 300–600 Da range.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phytosterols(smiles: str):
    """
    Determines if a molecule is a phytosterol based on its SMILES string.
    
    The algorithm:
      1. Parse the SMILES string.
      2. Extract the ring information and filter for rings that are made entirely of carbon.
      3. Build a connectivity graph among these rings (two rings are fused if they share >=2 atoms).
      4. For each connected component (candidate steroid nucleus), require:
             - Exactly four rings.
             - Among these four, one must be a 5-membered ring and three must be 6-membered.
             - The union of atoms in the nucleus (counting only carbons) is ~17 atoms.
      5. Check that at least one nucleus atom has a non‐nucleus neighbor that is part of a chain
         having at least 4 consecutive carbon atoms (candidate side chain).
      6. Check that the molecule contains at least one hydroxyl group (–OH).
      7. Check that the overall molecular weight is between 300 and 600 Da.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a phytosterol, False otherwise.
        str: Explanation (reason) for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if not rings:
        return False, "No rings found; not a steroid structure"
    
    # Filter rings: only consider rings made entirely of carbon atoms.
    pure_carbon_rings = []
    for ring in rings:
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            pure_carbon_rings.append(ring)
    if not pure_carbon_rings:
        return False, "No pure carbon rings found; steroid nucleus unlikely"
    
    # Build a graph among pure carbon rings: nodes are ring indices;
    # two rings are connected (fused) if they share at least 2 atoms.
    n_rings = len(pure_carbon_rings)
    ring_graph = {i: [] for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(pure_carbon_rings[i])
        for j in range(i+1, n_rings):
            set_j = set(pure_carbon_rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_graph[i].append(j)
                ring_graph[j].append(i)
                
    # Helper: DFS to get connected component in the ring graph
    def dfs(start, visited):
        stack = [start]
        component = set()
        while stack:
            node = stack.pop()
            if node not in visited:
                visited.add(node)
                component.add(node)
                for neigh in ring_graph[node]:
                    if neigh not in visited:
                        stack.append(neigh)
        return component

    visited = set()
    nucleus_found = False
    nucleus_atom_indices = set()
    nucleus_reason = ""
    
    # Look in each connected component for a candidate steroid nucleus
    for i in range(n_rings):
        if i in visited:
            continue
        comp = dfs(i, visited)
        # For a classic steroid nucleus, we expect exactly 4 fused rings.
        if len(comp) != 4:
            continue
        # Count ring sizes among these 4 rings (should be one 5-membered and three 6-membered).
        count5 = 0
        count6 = 0
        comp_atoms = set()
        for idx in comp:
            ring_atoms = pure_carbon_rings[idx]
            comp_atoms.update(ring_atoms)
            if len(ring_atoms) == 5:
                count5 += 1
            elif len(ring_atoms) == 6:
                count6 += 1
        if count5 == 1 and count6 == 3:
            # Check that the fused nucleus has ~17 carbon atoms.
            nC = sum(1 for idx in comp_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
            if 15 <= nC <= 19:
                nucleus_found = True
                nucleus_atom_indices = comp_atoms
                nucleus_reason = (f"Fused steroid nucleus: {len(comp)} rings (1 five-membered, 3 six-membered) "
                                   f"with {nC} carbon atoms.")
                break
    if not nucleus_found:
        return False, "Molecule does not contain a characteristic steroid nucleus (4 fused rings with required sizes)"
    
    # Check for an aliphatic side chain.
    # We assume that at least one atom in the nucleus has a neighbor that is not in the nucleus.
    # We then follow that neighbor (only counting carbons not already in the nucleus) to see if the chain length >= 4.
    def chain_length(start_idx, nucleus_set):
        # Simple DFS that only follows consecutive sp3 carbons not in the nucleus.
        seen = set()
        max_length = 0
        stack = [(start_idx, 1)]
        while stack:
            current, length = stack.pop()
            max_length = max(max_length, length)
            seen.add(current)
            for nbr in mol.GetAtomWithIdx(current).GetNeighbors():
                nbr_idx = nbr.GetIdx()
                # Only consider carbon atoms, not already in the nucleus.
                if nbr.GetAtomicNum() != 6 or nbr_idx in nucleus_set:
                    continue
                # Only follow sp3 (for an aliphatic chain) and avoid aromatic rings.
                if nbr.GetHybridization().name != "SP3":
                    continue
                if nbr_idx not in seen:
                    stack.append((nbr_idx, length+1))
        return max_length

    side_chain_found = False
    for idx in nucleus_atom_indices:
        atom = mol.GetAtomWithIdx(idx)
        # Check neighbors not in nucleus and that are carbon
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx not in nucleus_atom_indices and nbr.GetAtomicNum() == 6:
                # See how long is the chain starting from this neighbor.
                if chain_length(nbr_idx, nucleus_atom_indices) >= 4:
                    side_chain_found = True
                    break
        if side_chain_found:
            break
    if not side_chain_found:
        return False, "Steroid nucleus found but no significant aliphatic side chain attached"

    # Check for at least one hydroxyl group (-OH)
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    if hydroxyl_pattern is None:
        return False, "Error creating hydroxyl SMARTS pattern"
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Molecule lacks a hydroxyl group typical of phytosterols"
    
    # Check molecular weight (typical range for sterols ~300-600 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 600:
        return False, f"Molecular weight ({mol_wt:.1f} Da) is outside the typical sterol range for phytosterols"
    
    return True, ("Molecule exhibits a fused steroid nucleus (" + nucleus_reason +
                   "), an aliphatic side chain of sufficient length, a hydroxyl group, and an appropriate molecular weight, "
                   "consistent with a phytosterol.")

# Example usage (for testing purposes)
if __name__ == "__main__":
    # Example: ergosta-5,7-dien-3beta-ol (should be classified as phytosterol)
    test_smiles = "[H][C@@]1(CC[C@@]2([H])C3=CC=C4C[C@@H](O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CC[C@H](C)C(C)C"
    result, reason = is_phytosterols(test_smiles)
    print(result, reason)