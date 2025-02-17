"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Steroid saponin
Definition: Any saponin derived from a hydroxysteroid.
A steroid saponin is defined here as a molecule that (a) contains a hydroxysteroid nucleus—that is,
a fused tetracyclic core (approximately 15–19 atoms in total) composed of four candidate rings (3 six‐membered and 1 five‐membered, with most atoms being carbon) 
in which at least one atom of the core bears a hydroxyl (–OH) substituent—and (b) contains at least one sugar moiety.
A sugar moiety is approximated by detecting at least one 5‑ or 6‑membered ring that contains ≥2 oxygen atoms.
"""

from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.

    Steps:
      1. Parse the SMILES into an RDKit Mol and add explicit hydrogens.
      2. Identify candidate rings for the steroid nucleus:
         - Look at all rings of size 5 or 6.
         - For a 5-membered ring, require that at least 4 atoms are carbon.
           For a 6-membered ring, require that at least 5 atoms are carbon.
      3. Build connectivity among candidate rings. Two rings are considered fused if they share at least 2 atoms.
      4. From the connected clusters, search for a fused component comprising at least 4 rings
         and whose union of atoms is between 15 and 19 in number.
      5. Check that at least one atom of the fused system has a hydroxyl (-OH) substituent.
      6. Look for at least one sugar-like ring: a 5- or 6-membered ring that contains at least 2 oxygen atoms.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        (bool, str): (True, reason) if the molecule meets the steroid saponin criteria,
                     (False, reason) otherwise.
    """
    # Step 1. Parse SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Get ring information (all rings as tuples of atom indices) from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    
    # Step 2. Identify candidate rings for the steroid core.
    candidate_rings = []
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Count carbon atoms in the ring.
        carbon_count = sum(1 for atom in atoms if atom.GetAtomicNum() == 6)
        # Relax condition: For a 5-membered ring require at least 4 carbons;
        # for a 6-membered ring require at least 5 carbons.
        if (len(ring) == 5 and carbon_count >= 4) or (len(ring) == 6 and carbon_count >= 5):
            candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 4:
        return False, "Not enough candidate rings to form steroid nucleus"
    
    # Step 3. Build connectivity among candidate rings.
    # Two rings are fused if they share at least 2 atoms.
    n = len(candidate_rings)
    adj = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                adj[i].append(j)
                adj[j].append(i)
    
    # Step 4. Use DFS to find connected components of candidate rings.
    def dfs(start, visited):
        stack = [start]
        comp = set()
        while stack:
            cur = stack.pop()
            if cur in visited:
                continue
            visited.add(cur)
            comp.add(cur)
            for nb in adj[cur]:
                if nb not in visited:
                    stack.append(nb)
        return comp
    
    visited = set()
    fused_core_atoms = None
    # Look for any connected component (cluster of fused candidate rings)
    # that has at least 4 rings and a union of atoms between 15 and 19.
    for i in range(n):
        if i not in visited:
            comp = dfs(i, visited)
            if len(comp) >= 4:
                union_atoms = set()
                for idx in comp:
                    union_atoms |= candidate_rings[idx]
                if 15 <= len(union_atoms) <= 19:
                    fused_core_atoms = union_atoms
                    break
    if fused_core_atoms is None:
        return False, "No fused tetracyclic system consistent with a steroid nucleus found"
    
    # Step 5. Check that at least one atom in the steroid core bears a hydroxyl (-OH) substituent.
    hydroxyl_found = False
    for atom_idx in fused_core_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        # Look at neighbors for an -OH group.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen
                # Check if this oxygen is bound to a hydrogen.
                if any(nbr2.GetAtomicNum() == 1 for nbr2 in nbr.GetNeighbors()):
                    hydroxyl_found = True
                    break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "Steroid nucleus found, but no hydroxyl substituent detected on it"
    
    # Step 6. Look for a sugar-like ring. For our purposes, a sugar ring is defined as a 5- or
    # 6-membered ring that contains at least 2 oxygen atoms.
    sugar_found = False
    for ring in all_rings:
        if len(ring) in (5,6):
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxygen_count >= 2:
                sugar_found = True
                break
    if not sugar_found:
        return False, "No sugar-like ring (glycoside) found"
    
    return True, "Molecule contains a fused hydroxysteroid nucleus with a sugar moiety (steroid saponin)"

# Example usage (uncomment for testing):
# test_smiles = "O(C1C(C2(C(C3C(C4(C(=CC3)CC(O)CC4OC5OC(C(OC6OC(C(O)C(O)C6O)CO)C(O)C5O)C)CC2)C1)C)C(C(O)CCC(C)C)C)C7OC(C(O)C(O)C7O)CO"
# result, reason = is_steroid_saponin(test_smiles)
# print(result, reason)