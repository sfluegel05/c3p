"""
Classifies: CHEBI:61655 steroid saponin
"""
#!/usr/bin/env python
"""
Classifies: Steroid saponin
Definition: Any saponin derived from a hydroxysteroid.
A steroid saponin is defined here as a molecule that (a) contains a hydroxysteroid nucleus—that is,
a fused tetracyclic (approximately 17 carbon) core formed by four fused rings (3 six‐membered and 1 five‐membered) 
in which at least one of the core carbons is hydroxylated—and (b) contains at least one sugar moiety.
Here the steroid nucleus is approximated by searching for candidate rings (5- or 6-membered rings composed entirely 
of carbon atoms) that are fused (i.e. share at least a full bond, here an intersection of at least 2 atoms). 
A connected set of four such rings with 15–19 unique atoms is taken to represent the steroid core.
The sugar moiety is approximated by detecting at least one 5- or 6-membered ring that contains ≥2 oxygen atoms.
"""

from rdkit import Chem

def is_steroid_saponin(smiles: str):
    """
    Determines if a molecule is a steroid saponin based on its SMILES string.
    
    Steps:
      1. Parse the SMILES into an RDKit Mol.
      2. Search for a fused tetracyclic core: look for candidate rings (5- or 6-membered and all carbons)
         and use a connectivity search (via DFS) to find a connected set of ≥4 rings that share bonds (≥2 atoms overlap).
         We then require that the total number of unique atoms in that fused system is 15–19 (approximating a steroid nucleus).
      3. Verify that at least one atom within the steroid core has a hydroxyl (-OH) substituent.
      4. Look for at least one sugar-like ring: a 5- or 6-membered ring containing at least 2 oxygen atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): (True, reason) if the molecule meets the criteria,
                     (False, reason) if it does not.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # tuples of atom indices

    # Step 1. Identify candidate rings:
    candidate_rings = []
    for ring in all_rings:
        if len(ring) in (5, 6):
            # Check if every atom in the ring is carbon (atomic number 6)
            if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
                candidate_rings.append(set(ring))
    
    if len(candidate_rings) < 4:
        return False, "Not enough candidate rings to form steroid nucleus"
    
    # Step 2. Build connectivity among candidate rings:
    # Two rings are fused if they share at least two atoms (i.e. a full bond)
    n = len(candidate_rings)
    adj = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(candidate_rings[i].intersection(candidate_rings[j])) >= 2:
                adj[i].append(j)
                adj[j].append(i)
    
    # DFS to obtain connected components (clusters of fused rings)
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
    # Look for any connected component with at least 4 rings whose union has between 15 and 19 atoms.
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
    
    # Step 3. Check that at least one atom of the steroid core has a hydroxyl (-OH) substituent.
    hydroxyl_found = False
    for atom_idx in fused_core_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:  # oxygen
                # Check if oxygen is attached to a hydrogen (i.e., -OH)
                for subnbr in nbr.GetNeighbors():
                    if subnbr.GetAtomicNum() == 1:
                        hydroxyl_found = True
                        break
            if hydroxyl_found:
                break
        if hydroxyl_found:
            break
    if not hydroxyl_found:
        return False, "Steroid nucleus found, but no hydroxyl substituent detected on it"
    
    # Step 4. Look for a sugar-like ring.
    # A sugar-like ring is a 5- or 6-membered ring that contains at least 2 oxygen atoms.
    sugar_found = False
    for ring in all_rings:
        if len(ring) in (5, 6):
            oxygen_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxygen_count >= 2:
                sugar_found = True
                break
    if not sugar_found:
        return False, "No sugar-like ring (glycoside) found"
    
    return True, "Molecule contains a fused hydroxysteroid nucleus with a sugar moiety (steroid saponin)"

# For example, you might test using:
# test_smiles = "O(C1C(C2(C(C3C(C4(C(=CC3)CC(O)CC4OC5OC(C(OC6OC(C(O)C(O)C6O)CO)C(O)C5O)C)CC2)C1)C)C(C(O)CCC(C)C)C)C7OC(C(O)C(O)C7O)CO"
# result, reason = is_steroid_saponin(test_smiles)
# print(result, reason)