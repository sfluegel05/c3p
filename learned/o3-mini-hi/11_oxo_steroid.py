"""
Classifies: CHEBI:47787 11-oxo steroid
"""
#!/usr/bin/env python
"""
Classifies: 11-oxo steroid
Definition: An oxo steroid is any steroid-like molecule that contains a ketone (C=O) at position 11.
This improved heuristic proceeds as follows:
  1. Find a fused ring system with four rings. Traditionally the steroid nucleus is three six-membered rings and one five-membered ring.
     To accommodate variants (e.g. D-homo steroids) we allow the non-six-membered ring to be either 5- or 7-membered, while the other three must be 6-membered.
  2. The four rings must be interconnected (fused) and the union of their atoms should fall in a characteristic range.
  3. Identify the “key” six-membered ring that is most fused with the non-six-membered (D or D-homo) ring. Then search that ring, or atoms adjacent to it,
     for a ketone group (C=O). Only if such a ketone is found do we call the molecule an 11-oxo steroid.
Note:
  This method is heuristic and will not capture every exotic steroid but aims to reduce false positives and negatives.
"""
from rdkit import Chem
from itertools import combinations

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string with an improved heuristic.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        bool: True if the molecule is classified as an 11-oxo steroid, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    ring_info = mol.GetRingInfo()
    # Allow rings of sizes 5, 6, or 7 for nucleus detection.
    candidate_rings = [set(r) for r in ring_info.AtomRings() if len(r) in (5,6,7)]
    if len(candidate_rings) < 4:
        return False, "Not enough 5-, 6-, or 7-membered rings for steroid nucleus."

    nucleus_atoms = None
    nucleus_ring_combo = None
    # Loop over sets of 4 rings.
    for ring_combo in combinations(candidate_rings, 4):
        # For a classical steroid nucleus, we expect three rings to be six-membered
        # and one ring to be the 'D-ring' which can be either 5- or 7-membered.
        count_non6 = sum(1 for r in ring_combo if len(r) != 6)
        count6 = sum(1 for r in ring_combo if len(r) == 6)
        if count6 != 3 or count_non6 != 1:
            continue

        # Check that every ring in the combination is fused to at least one other.
        nodes = list(ring_combo)
        adj = {i: set() for i in range(4)}
        for i in range(4):
            for j in range(i+1,4):
                if nodes[i] & nodes[j]:
                    adj[i].add(j)
                    adj[j].add(i)
        visited = set()
        stack = [0]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            stack.extend(adj[node] - visited)
        if len(visited) != 4:
            continue  # Rings are not fully fused.

        # Compute the union of atom indices in these rings.
        union_atoms = set()
        for r in ring_combo:
            union_atoms.update(r)
        # Require the union to be in a narrower range (e.g. 15 to 25 atoms).
        if not (15 <= len(union_atoms) <= 25):
            continue
        # Require that the nucleus has a sufficient number of carbons.
        n_carbons = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbons < 12:
            continue

        nucleus_atoms = union_atoms
        nucleus_ring_combo = ring_combo
        break

    if nucleus_atoms is None:
        return False, "No fused steroid nucleus (4 rings with appropriate sizes and connectivity) found."

    # Identify the unique non-six-membered ring (the candidate D-ring)
    D_ring = None
    six_membered = []
    for r in nucleus_ring_combo:
        if len(r) != 6:
            D_ring = r
        else:
            six_membered.append(r)
    if D_ring is None or not six_membered:
        # Fallback: if our expected pattern did not arise, use the whole nucleus.
        candidate_key_ring = nucleus_atoms
    else:
        # Among the six-membered rings, select the one sharing most atoms with the D-ring.
        candidate_key_ring = None
        max_shared = 0
        for r in six_membered:
            shared = len(r & D_ring)
            if shared > max_shared:
                max_shared = shared
                candidate_key_ring = r
        # If none of the six-membered rings touches the D ring, fallback to the entire nucleus.
        if candidate_key_ring is None:
            candidate_key_ring = nucleus_atoms

    # Now search for a ketone (C=O) group.
    # Instead of using a broad SMARTS, we manually check ring atoms for a direct double-bond to O.
    # Iterate over atoms in candidate_key_ring.
    for idx in candidate_key_ring:
        atom = mol.GetAtomWithIdx(idx)
        # Check is it a carbon?
        if atom.GetAtomicNum() != 6:
            continue
        # Look at bonds: if any double bond to an oxygen is found, we have a ketone on the ring.
        for bond in atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2:  # double bond
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    return True, ("Steroid nucleus detected and a ketone (C=O) directly attached to the key ring "
                                  "(consistent with an 11-oxo steroid).")
        # Also check: even if the carbon itself is not in a carbonyl, it might be adjacent to one that is.
        # (This allows for minor variations of connectivity.)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            for bond in nbr.GetBonds():
                # If the neighbor has a double-bond to oxygen and is directly bonded to an atom in the key ring
                if bond.GetBondTypeAsDouble() == 2:
                    other = bond.GetOtherAtom(nbr)
                    if other.GetAtomicNum() == 8:
                        return True, ("Steroid nucleus detected and a ketone (C=O) on a neighbor of the key ring "
                                      "(consistent with an 11-oxo steroid).")
    
    # If no ketone was found on or adjacent to the candidate ring, then it is not an 11-oxo steroid.
    return False, ("Ketone group not found on or immediately attached to the key ring of the steroid nucleus "
                   "(likely not a true 11-oxo steroid).")

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        # True positives (from provided examples)
        "O=C1C2=C([C@@]3(C(=O)C[C@@H]([C@]3(C1)C)C(=C)CCC(=O)O)C)[C@@H](O)C[C@@H]4[C@@]2(CCC(C4(C)C)=O)",  # 20(21)-dehydrolucidenic acid A
        "C1[C@@]2([C@@]([C@@]3(C(C1)=CC(CC3)=O)C)(C(C[C@@]4(C)[C@@]2([H])CC[C@@H]4O)=O)[H])[H]",  # 11-oxotestosterone
        "[H][C@]12CC[C@@]3([H])[C@]4([H])CCC\\C(=C/CO)[C@@]4(C)CC(=O)[C@]3([H])[C@@]1(C)CC[C@@H](O)C2",  # example previously missed (should be TP)
        "OC(=O)C(C)C(=C)CC[C@@H](C)[C@H]1CC[C@H]2C3=C(C(=O)C[C@]12C)[C@@]1(C)CCC(=O)[C@@H](C)[C@@H]1CC3=O",  # zhankuic acid A methyl ester (TP)
        # False positive example (should not be classified as 11-oxo steroid)
        "O=C1C2=C(O)C=3C(=O)C=C(OC)C(C3C=C2C(=O)C4=C1C5=C(C(=O)OC5)C(=C4)OC)=O",  # Paramagnetoquinone B
    ]
    for smi in test_smiles:
        result, reason = is_11_oxo_steroid(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-----")