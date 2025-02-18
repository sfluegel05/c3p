"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-Delta(4) steroid
Definition: A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position,
i.e. a steroid (with a fused tetracyclic nucleus) that contains an α,β-unsaturated ketone 
(enone) in one of its six-membered rings.
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    The classification requires:
      (a) A fused steroid nucleus – defined as a connected system of rings (only considering rings of size 5 or 6)
          containing at least 4 rings (typically three six-membered and one five-membered).
      (b) Within one of the six-membered rings of that fused system, an enone motif is present; that is,
          a carbon (in the ring) that has a double bond to an oxygen (C=O), and via a single bond is connected 
          to a second carbon in the ring that, in turn, participates in a C=C double bond with a third carbon (within the same ring).
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # ---------------------------
    # 1. Identify candidate rings (only 5- or 6-membered rings).
    # ---------------------------
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    candidate_rings = []
    for ring in all_rings:
        if len(ring) in (5, 6):
            candidate_rings.append(set(ring))
    if not candidate_rings:
        return False, "No candidate 5- or 6-membered rings found for steroid nucleus"

    # ---------------------------
    # 2. Build a fused ring graph.
    # Two rings will be considered "fused" if they share at least 2 atoms.
    # ---------------------------
    num = len(candidate_rings)
    graph = {i: set() for i in range(num)}
    for i in range(num):
        for j in range(i+1, num):
            if len(candidate_rings[i] & candidate_rings[j]) >= 2:
                graph[i].add(j)
                graph[j].add(i)

    # Find connected components in the ring graph.
    visited = set()
    fused_components = []
    for i in range(num):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                current = stack.pop()
                if current in comp:
                    continue
                comp.add(current)
                for nbr in graph[current]:
                    if nbr not in comp:
                        stack.append(nbr)
            visited |= comp
            fused_components.append(comp)

    # ---------------------------
    # 3. Look for a fused tetracyclic nucleus.
    # We require at least 4 rings in one connected component with at least three six‐membered rings and one five‐membered ring.
    # ---------------------------
    steroid_component = None
    for comp in fused_components:
        comp_rings = [candidate_rings[i] for i in comp]
        count6 = sum(1 for ring in comp_rings if len(ring) == 6)
        count5 = sum(1 for ring in comp_rings if len(ring) == 5)
        if len(comp_rings) >= 4 and count6 >= 3 and count5 >= 1:
            # Found a candidate steroid nucleus
            steroid_component = comp_rings  # list of sets of atom indices
            break
    if steroid_component is None:
        return False, "Steroid nucleus not found (insufficient fused ring system with 3 six-membered and 1 five-membered rings)"

    # Build a union of all atom indices in the fused steroid nucleus:
    steroid_atoms = set()
    for ring in steroid_component:
        steroid_atoms.update(ring)

    # ---------------------------
    # 4. Look for an enone functionality (α,β-unsaturated ketone) in one six-membered ring of the nucleus.
    # ---------------------------
    # We restrict our search to six-membered rings that are part of the steroid nucleus.
    six_membered_rings = [ring for ring in steroid_component if len(ring)==6]
    enone_found = False
    # Iterate through each six-membered ring
    for ring in six_membered_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Consider only carbon atoms
            if atom.GetAtomicNum() != 6:
                continue
            # Look for a carbonyl: check if this carbon is double-bonded to at least one oxygen.
            carbonyl_oxygens = []
            for bond in atom.GetBonds():
                # check if bond is double and neighbor is oxygen
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    nbr = bond.GetOtherAtom(atom)
                    if nbr.GetAtomicNum() == 8:
                        carbonyl_oxygens.append(nbr.GetIdx())
            if not carbonyl_oxygens:
                continue  # no carbonyl at this carbon

            # Now search for a beta carbon: a neighbor (in the same ring) via a single bond.
            for bond in atom.GetBonds():
                nbr = bond.GetOtherAtom(atom)
                if nbr.GetIdx() not in ring:
                    continue  # must be part of the same ring
                if nbr.GetAtomicNum() != 6:
                    continue  # ensure neighbor is carbon
                # Avoid the carbonyl oxygen connection
                if nbr.GetIdx() in carbonyl_oxygens:
                    continue

                # Now check if this neighbor (the beta carbon) participates in a C=C double bond to another carbon in the same ring.
                for bond2 in nbr.GetBonds():
                    if bond2.GetBondType() != Chem.rdchem.BondType.DOUBLE:
                        continue
                    gamma = bond2.GetOtherAtom(nbr)
                    if gamma.GetIdx() == atom.GetIdx():
                        continue  # skip if going back to original carbonyl carbon
                    if gamma.GetIdx() not in ring:
                        continue  # ensure gamma carbon is in the same ring
                    if gamma.GetAtomicNum() != 6:
                        continue
                    # Found a candidate enone: carbonyl carbon (with double bond O) connected (single bond) to a beta carbon that participates in a double bond with gamma.
                    enone_found = True
                    break
                if enone_found:
                    break
            if enone_found:
                break
        if enone_found:
            break

    if not enone_found:
        return False, "Enone motif not found in a six-membered ring within the steroid nucleus"
    
    return True, "Contains a fused steroid nucleus (3 six-membered and 1 five-membered fused rings) with a 3-oxo/Δ(4) enone motif."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with (20S)-20-hydroxypregn-4-en-3-one
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CCC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)O"
    result, reason = is_3_oxo_Delta_4__steroid(test_smiles)
    print(result, reason)