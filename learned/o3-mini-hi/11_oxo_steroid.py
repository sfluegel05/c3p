"""
Classifies: CHEBI:47787 11-oxo steroid
"""
#!/usr/bin/env python
"""
Classifies: 11-oxo steroid
Definition: Any oxo steroid that has an oxo substituent (ketone) at position 11.
This heuristic first detects a steroid nucleus – defined as four fused rings 
(three six-membered and one five-membered) whose union of atoms falls (relaxed) into
a count range—and then identifies the “key” six-membered ring (that is fused to the 
five-membered ring, corresponding to the C ring where C11 is located). Only a ketone 
(C=O) found on or immediately attached to that key ring is accepted.
Note:
  This method is heuristic and may still miss very exotic cases.
"""

from rdkit import Chem
from itertools import combinations

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    The approach is to first search for a steroid nucleus, defined as a set of 4 fused rings 
    (exactly one five-membered and three six-membered rings) whose union of atom indices falls 
    in a relaxed range. Next, among the nucleus rings the five-membered ring (the D ring) is identified;
    then among the six-membered rings, the one sharing the most atoms with the five-membered ring 
    is selected as the candidate C ring (where the 11-oxo is expected). Finally, we look for a ketone 
    (C=O) group whose carbonyl carbon is either in that candidate ring or directly adjacent to it.
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if classified as 11-oxo steroid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Get ring information and filter rings to sizes 5 or 6.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    candidate_rings = [set(r) for r in all_rings if len(r) in (5, 6)]
    if len(candidate_rings) < 4:
        return False, "Not enough 5- or 6-membered rings to form a steroid nucleus."
    
    nucleus_atoms = None
    nucleus_ring_combo = None
    # We relax union count to between 14 and 24 atoms and also require at least 12 carbon atoms.
    for ring_combo in combinations(candidate_rings, 4):
        count5 = sum(1 for r in ring_combo if len(r) == 5)
        count6 = sum(1 for r in ring_combo if len(r) == 6)
        if count5 != 1 or count6 != 3:
            continue  # Must have exactly one 5-membered and three 6-membered rings.
        # Check connectivity among these rings: build a graph where each ring is a node.
        nodes = list(ring_combo)
        adj = {i: set() for i in range(4)}
        for i in range(4):
            for j in range(i+1, 4):
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
            continue  # Rings not all fused.
        # Compute union of atoms and check total count.
        union_atoms = set()
        for r in ring_combo:
            union_atoms.update(r)
        if not (14 <= len(union_atoms) <= 24):
            continue
        # Count carbons in the union.
        n_carbons = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbons < 12:
            continue
        # We have a candidate nucleus.
        nucleus_atoms = union_atoms
        nucleus_ring_combo = ring_combo
        break
    if nucleus_atoms is None:
        return False, "No steroid nucleus (4 fused rings with appropriate sizes) found."
    
    # Identify the unique five-membered ring among the nucleus rings.
    five_membered = None
    six_membered = []
    for r in nucleus_ring_combo:
        if len(r) == 5:
            five_membered = r
        else:
            six_membered.append(r)
    if five_membered is None or not six_membered:
        # fallback if not identified as expected.
        candidate_key_ring = None
    else:
        # Among the six-membered rings, choose the one sharing the most atoms with the five-membered ring.
        candidate_key_ring = None
        max_shared = 0
        for r in six_membered:
            shared = len(r & five_membered)
            if shared > max_shared:
                max_shared = shared
                candidate_key_ring = r
        # It is possible that no six-membered ring is fused with the five-membered ring;
        # in that case, we fallback to checking all nucleus atoms.
    
    # Prepare a SMARTS for a ketone group.
    # Here the pattern [CX3](=O) finds a carbonyl group (ketone) and we allow any connectivity.
    ketone_smarts = "[CX3](=O)"
    ketone_pat = Chem.MolFromSmarts(ketone_smarts)
    if ketone_pat is None:
        return None, None  # SMARTS creation failed.
    
    ketone_matches = mol.GetSubstructMatches(ketone_pat)
    if not ketone_matches:
        return False, "No ketone (C=O) group detected (required for 11-oxo steroid)."
    
    # Now check for a ketone that is ideally located on the candidate key ring
    # (i.e. the six-membered ring fused with the five-membered ring) or directly attached to it.
    for match in ketone_matches:
        # In our ketone pattern, match[0] is the carbonyl carbon.
        ketone_carbon = match[0]
        in_candidate = False
        if candidate_key_ring is not None:
            if ketone_carbon in candidate_key_ring:
                in_candidate = True
            else:
                # Check if any neighbor of the ketone carbon is in the candidate ring.
                atom = mol.GetAtomWithIdx(ketone_carbon)
                if any(neighbor.GetIdx() in candidate_key_ring for neighbor in atom.GetNeighbors()):
                    in_candidate = True
        # If no candidate key ring was found, fallback to checking whether the ketone is on the nucleus 
        # or attached to any nucleus atom.
        else:
            if ketone_carbon in nucleus_atoms:
                in_candidate = True
            else:
                atom = mol.GetAtomWithIdx(ketone_carbon)
                if any(neighbor.GetIdx() in nucleus_atoms for neighbor in atom.GetNeighbors()):
                    in_candidate = True
        
        if in_candidate:
            return True, ("Steroid nucleus detected and ketone group found on or immediately attached "
                          "to the key nucleus ring (consistent with an 11-oxo steroid).")
    
    return False, ("Ketone group not located on or adjacent to the steroid nucleus's key ring "
                   "(likely not a true 11-oxo steroid).")

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "O=C1C2=C([C@@]3(C(=O)C[C@@H]([C@]3(C1)C)C(=C)CCC(=O)O)C)[C@@H](O)C[C@@H]4[C@@]2(CCC(C4(C)C)=O)",  # 20(21)-dehydrolucidenic acid A (TP)
        "C1[C@@]2([C@@]([C@@]3(C(C1)=CC(CC3)=O)C)(C(C[C@@]4(C)[C@@]2([H])CC[C@@H]4O)=O)[H])[H]",  # 11-oxotestosterone (TP)
        "[H][C@]12CC[C@@]3([H])[C@]4([H])CCC\\C(=C/CO)[C@@]4(C)CC(=O)[C@]3([H])[C@@]1(C)CC[C@@H](O)C2",  # Previously false negative (should be TP)
        "[H][C@]12C[C@@]1([H])[C@]1(C)[C@@H](O)C(=O)\\C(=C(\\C)C(=O)OC)[C@@]3([H])C1=C2C[C@@]1([H])[C@@]2(C)[C@]4([H])C[C@]4([H])[C@](O)(C\\C(C)=C\\C(=O)OC)[C@@]2([H])CC2=C(CO)C(=O)O[C@]312",  # cortistatin D (false positive)
    ]
    for smi in test_smiles:
        res, reason = is_11_oxo_steroid(smi)
        print("SMILES:", smi)
        print("Result:", res)
        print("Reason:", reason)
        print("-----")