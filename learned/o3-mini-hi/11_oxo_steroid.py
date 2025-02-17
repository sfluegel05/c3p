"""
Classifies: CHEBI:47787 11-oxo steroid
"""
#!/usr/bin/env python
"""
Classifies: 11-oxo steroid
Definition: A steroid is classified as an 11-oxo steroid if a fused ring nucleus (three six‐membered rings and one non‐six‐membered ring)
is detected and if a ketone (C=O) group is found either on or immediately adjacent to the fused nucleus.
This improved heuristic relaxes size constraints on the nucleus and broadens the ketone search region.
"""

from rdkit import Chem
from itertools import combinations

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.
    This function looks for a fused 4-ring steroid nucleus and then examines the nucleus
    and its immediate neighbors for a carbonyl group (C=O).
    
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
    # We are interested in rings of size 5, 6, or 7 (allowing for D-ring variants)
    candidate_rings = [set(r) for r in ring_info.AtomRings() if len(r) in (5, 6, 7)]
    if len(candidate_rings) < 4:
        return False, "Not enough 5-, 6-, or 7-membered rings for a steroid nucleus."

    nucleus_atoms = None
    nucleus_ring_combo = None

    # Loop over all candidate combinations of 4 rings
    # We require that 3 rings are 6-membered and 1 ring is non-6 (either 5- or 7-membered)
    for ring_combo in combinations(candidate_rings, 4):
        count_six = sum(1 for r in ring_combo if len(r) == 6)
        count_non6 = 4 - count_six
        if count_six != 3 or count_non6 != 1:
            continue

        # Check that the rings are fused (each ring must share at least one atom with at least one other)
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
            continue  # not fully fused

        # Compute the union of the atoms in the four rings.
        union_atoms = set()
        for r in ring_combo:
            union_atoms.update(r)
        # Instead of a strict range of 15-25 atoms (which might miss borderline cases),
        # we relax the requirement a little.
        if not (14 <= len(union_atoms) <= 30):
            continue

        # Require that the nucleus has a minimum number of carbon atoms.
        n_carbons = sum(1 for idx in union_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if n_carbons < 12:
            continue

        # Accept this combination as the steroid nucleus.
        nucleus_atoms = union_atoms
        nucleus_ring_combo = ring_combo
        break

    if nucleus_atoms is None:
        return False, "No fused steroid nucleus (4 rings with appropriate sizes and connectivity) found."

    # Instead of focusing on one key ring, we expand the search region to include:
    # (a) all atoms in the nucleus and (b) all immediate neighbors of those atoms.
    search_region = set(nucleus_atoms)
    for idx in nucleus_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            search_region.add(nbr.GetIdx())

    # Now search for a ketone (C=O) group in the search region.
    # We look for a carbon atom double-bonded to an oxygen.
    # To be conservative, we check all atoms in the search region.
    for idx in search_region:
        atom = mol.GetAtomWithIdx(idx)
        if atom.GetAtomicNum() != 6:
            continue
        for bond in atom.GetBonds():
            # Check for a double bond (order 2) to an oxygen
            if bond.GetBondTypeAsDouble() == 2:
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetAtomicNum() == 8:
                    # We found a carbonyl within or immediately adjacent to the nucleus.
                    # Provide a message depending on whether the carbon is in the nucleus.
                    if idx in nucleus_atoms:
                        detail = ("Steroid nucleus detected and a ketone (C=O) directly attached "
                                  "to the nucleus (consistent with an 11-oxo steroid).")
                    else:
                        detail = ("Steroid nucleus detected and a ketone (C=O) immediately adjacent "
                                  "to the nucleus (consistent with an 11-oxo steroid).")
                    return True, detail

    return False, ("Ketone group not found on or adjacent to the fused steroid nucleus "
                   "(likely not a true 11-oxo steroid).")


# Example usage (testing some provided examples):
if __name__ == "__main__":
    test_smiles = [
        # True positives
        "O=C1C2=C([C@@]3(C(=O)C[C@@H]([C@]3(C1)C)C(=C)CCC(=O)O)C)[C@@H](O)C[C@@H]4[C@@]2(CCC(C4(C)C)=O)",  # 20(21)-dehydrolucidenic acid A
        "O=C1C2=C([C@@]3(C(=O)C[C@@H]([C@]3(C1)C)C(=C)CCC(=O)OC)C)[C@@H](O)C[C@@H]4[C@@]2(CCC(C4(C)C)=O)",  # methyl 20(21)-dehydrolucidenate A
        "C1[C@@]2([C@@]([C@@]3(C(C1)=CC(CC3)=O)C)(C(C[C@@]4(C)[C@@]2([H])CC[C@@H]4O)=O)[H])[H]",  # 11-oxotestosterone
        "C[C@]12CC(=O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@]2(O)C(=O)C=O",  # 21-dehydrocortisone
        # False positive sample (should not be classified as 11-oxo steroid)
        "O=C1C2=C(O)C=3C(=O)C=C(OC)C(C3C=C2C(=O)C4=C1C5=C(C(=O)OC5)C(=C4)OC)=O",  # Paramagnetoquinone B
        # False negative example (previously missed)
        "[H][C@]12CC[C@@]3([H])[C@]4([H])CCC\\C(=C/CO)[C@@]4(C)CC(=O)[C@]3([H])[C@@]1(C)CC[C@@H](O)C2",  # 3alpha,21-dihydroxy-D-homo-5beta-pregn-17a(20)-en-11-one
    ]

    for smi in test_smiles:
        result, reason = is_11_oxo_steroid(smi)
        print("SMILES:", smi)
        print("Result:", result)
        print("Reason:", reason)
        print("-----")