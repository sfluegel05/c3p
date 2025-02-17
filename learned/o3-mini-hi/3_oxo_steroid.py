"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition: Any oxo steroid where an oxo substituent (a ketone group, i.e. C=O) is located on a six‐membered ring of the steroid nucleus.
Heuristic used in this implementation:
  1. Identify all rings of size 5 or 6.
  2. Build a graph of rings that are fused (sharing at least 2 atoms) and extract the largest connected component.
  3. Consider this fused system as the steroid “nucleus” if it contains at least 4 rings and at least 3 are six‐membered.
  4. Additionally, require that the number of (carbon) atoms in the nucleus is in the typical range (15 to 21).
  5. Finally, search for a ketone group (a double-bonded oxygen) on one of the six‐membered rings in the nucleus.
If any of these checks fail, the molecule is rejected.
"""

from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is (heuristically) defined as a molecule that contains a fused ring (steroid) nucleus
    with at least 4 rings (of which at least 3 are six‐membered) and the total number of spanned carbons is in a typical range.
    In this nucleus, at least one six‐membered ring has a ketone (C=O) substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-oxo steroid, False otherwise.
        str : A reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information (each ring as a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found; not a steroid"
    
    # Consider only rings of size 5 or 6.
    candidate_rings = [ring for ring in ring_info if len(ring) in (5, 6)]
    if len(candidate_rings) < 4:
        return False, "Less than 4 candidate rings (size 5 or 6); not a typical steroid nucleus"

    # Build a graph of which candidate rings are fusion-connected (sharing at least 2 atoms)
    n = len(candidate_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i+1, n):
            if len(set(candidate_rings[i]).intersection(candidate_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components in the ring graph
    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                current = stack.pop()
                if current not in comp:
                    comp.add(current)
                    visited.add(current)
                    stack.extend(ring_graph[current] - comp)
            components.append(comp)
    
    # Use the largest connected component as our nucleus candidate.
    largest_comp = max(components, key=lambda c: len(c))
    if len(largest_comp) < 4:
        return False, "Largest fused ring system does not have at least 4 rings; not a steroid nucleus"
    
    # Count how many of the rings in the nucleus are six-membered and five-membered.
    six_count = 0
    five_count = 0
    nucleus_atoms = set()
    for idx in largest_comp:
        ring = candidate_rings[idx]
        nucleus_atoms.update(ring)
        if len(ring) == 6:
            six_count += 1
        elif len(ring) == 5:
            five_count += 1

    if six_count < 3:
        return False, f"Fused ring system found but only {six_count} six-membered ring(s) detected, expected at least 3"
    # Additionally, check if the total count of carbon atoms (from the nucleus atoms) is within a reasonable range.
    nucleus_carbon_count = sum(1 for atom in mol.GetAtoms() 
                               if atom.GetAtomicNum() == 6 and atom.GetIdx() in nucleus_atoms)
    if nucleus_carbon_count < 15 or nucleus_carbon_count > 21:
        return False, f"Fused ring nucleus has {nucleus_carbon_count} carbons, which is out of range for a typical steroid nucleus"

    # Now, search for a ketone group (C=O) on one of the six-membered rings of the nucleus.
    # We look at each candidate ring (of size 6) in the nucleus.
    for idx in largest_comp:
        ring = candidate_rings[idx]
        if len(ring) == 6:
            # For each carbon in this six-membered ring
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    # Check bonds for a double bond to oxygen.
                    for bond in atom.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                # Also check that this carbon is part of the nucleus.
                                return True, ("Found ketone group on a six-membered ring of the steroid nucleus "
                                              "(likely at position 3)")
    return False, ("No ketone group found on any six-membered ring within the fused steroid nucleus")

# Example usage (testing):
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO",  # methyl prednisolone-16alpha-carboxylate
        "C[C@]12CC[C@H]3[C@@H](CC=C4CC(=O)CC[C@]34C)[C@@H]1CCC2=O",  # androst-5-ene-3,17-dione
        # False negative example (previously missed due to extra rings):
        "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])CC[C@@]1(C)[C@@]2([H])C[C@H]2OC(C)(C)O[C@@]12C(C)=O",  # algestone acetonide
        # A false positive candidate (should be rejected by the nucleus filters):
        "C[C@H](CCC(=O)NCCCC(O)=O)[C@H]1CC[C@H]2[C@H]3[C@H](CC[C@]12C)[C@@]4(C)CC[C@@H](O)C[C@H]4CC3=O"  # 7-oxo-gamma-aminoisobutyrolithocholic acid
    ]
    for sm in test_smiles:
        result, reason = is_3_oxo_steroid(sm)
        print(sm[:70], "->", result, "|", reason)