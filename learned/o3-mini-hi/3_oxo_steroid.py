"""
Classifies: CHEBI:47788 3-oxo steroid
"""
"""
Classifies: 3-oxo steroid
Definition: Any oxo steroid where a ketone (C=O) substituent is located on a six‐membered ring of the steroid nucleus.
Heuristics used in this implementation:
  1. Identify all rings of size 5 or 6.
  2. Build a graph of fused rings (rings that share at least 2 atoms) and choose the largest connected component as the candidate steroid nucleus.
  3. Verify that the candidate nucleus contains at least 4 rings, that it contains at least one five–membered ring (typical of steroid D-ring), and that the overall number of carbon atoms in the nucleus is in a typical range.
  4. Relax the requirement on the number of six–membered rings if the nucleus is small (i.e. a “nor–steroid”).
  5. Search for a ketone (C=O) group on any six–membered ring in the nucleus.
If any check fails the molecule is rejected.
"""

from rdkit import Chem

def is_3_oxo_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo steroid based on its SMILES string.
    A 3-oxo steroid is (heuristically) defined as a molecule that contains a fused steroid nucleus
    with 4 or more rings (with at least one being five–membered and at least 2–3 being six–membered, 
    depending on the overall carbon count) and that at least one six–membered ring bears a ketone (C=O) group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3-oxo steroid, False otherwise.
        str : Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get all rings (as tuples of atom indices) from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found; not a steroid"

    # Consider only rings of size 5 or 6.
    candidate_rings = [ring for ring in ring_info if len(ring) in (5, 6)]
    if len(candidate_rings) < 4:
        return False, "Less than 4 candidate rings (size 5 or 6); not a typical steroid nucleus"

    # Build a graph where nodes are candidate rings, and edges exist if they share at least 2 atoms.
    n = len(candidate_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if len(set(candidate_rings[i]).intersection(candidate_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)

    # Find connected components (each represents a fused ring system).
    visited = set()
    components = []
    for i in range(n):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                cur = stack.pop()
                if cur not in comp:
                    comp.add(cur)
                    visited.add(cur)
                    stack.extend(ring_graph[cur] - comp)
            components.append(comp)

    # Use the largest connected component as the candidate nucleus.
    largest_comp = max(components, key=lambda c: len(c))
    if len(largest_comp) < 4:
        return False, "Fused ring system does not have at least 4 rings; not a steroid nucleus"

    # Count six–membered and five–membered rings in the nucleus and collect all nucleus atom indices.
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

    # Count carbon atoms in the nucleus.
    nucleus_carbon_count = sum(1 for atom in mol.GetAtoms() 
                               if atom.GetAtomicNum() == 6 and atom.GetIdx() in nucleus_atoms)
    # Typical steroid nuclei (or nor-steroids) tend to have between 15 and 24 carbons.
    if nucleus_carbon_count < 15 or nucleus_carbon_count > 24:
        return False, f"Fused ring nucleus has {nucleus_carbon_count} carbons, which is out of the typical range (15-24)"

    # For a classical steroid, expect three six-membered rings and one five-membered ring.
    # Relax the requirement to 2 six-membered rings if the nucleus has fewer carbons (e.g. nor-steroids).
    required_six = 2 if nucleus_carbon_count <= 17 else 3
    if six_count < required_six:
        return False, f"Fused ring system has only {six_count} six-membered ring(s), expected at least {required_six}"
    if five_count < 1:
        return False, "Fused ring nucleus does not contain a five-membered ring typical of steroids"

    # Now, search for a ketone (C=O) on any six-membered ring within the nucleus.
    for idx in largest_comp:
        ring = candidate_rings[idx]
        if len(ring) == 6:
            for atom_idx in ring:
                atom = mol.GetAtomWithIdx(atom_idx)
                if atom.GetAtomicNum() == 6:
                    for bond in atom.GetBonds():
                        # Look for a double bond to an oxygen.
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            nbr = bond.GetOtherAtom(atom)
                            if nbr.GetAtomicNum() == 8:
                                return True, ("Found ketone group on a six-membered ring of the steroid nucleus "
                                              "(likely at position 3)")
    return False, "No ketone group found on any six-membered ring within the fused steroid nucleus"

# Example usage for testing:
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "[H][C@@]12CCC3=CC(=O)C=C[C@]3(C)[C@@]1([H])[C@@H](O)C[C@@]1(C)[C@@]2([H])C[C@@H](C(=O)OC)[C@]1(O)C(=O)CO",  # methyl prednisolone-16alpha-carboxylate
        "C[C@]12CC[C@H]3[C@@H](CC=C4CC(=O)CC[C@]34C)[C@@H]1CCC2=O",  # androst-5-ene-3,17-dione
        # Example that previously triggered false negative:
        "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)CCCC(C)C",  # B-norcholest-4-en-3-one
        # False positive candidate (should be rejected):
        "C[C@H](CCC(=O)NCCCC(O)=O)[C@H]1CC[C@H]2[C@H]3[C@H](CC[C@]12C)[C@@]4(C)CC[C@@H](O)C[C@H]4CC3=O"  # 7-oxo-gamma-aminoisobutyrolithocholic acid
    ]
    for sm in test_smiles:
        result, reason = is_3_oxo_steroid(sm)
        print(sm[:70], "->", result, "|", reason)