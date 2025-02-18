"""
Classifies: CHEBI:47909 3-oxo-Delta(4) steroid
"""
#!/usr/bin/env python3
"""
Classifies: 3-oxo-Delta(4) steroid
Definition: A 3-oxo steroid conjugated to a C=C double bond at the alpha,beta position,
i.e. a steroid (with a fused tetracyclic nucleus) that contains an α,β-unsaturated ketone 
(enone) moiety in one of its six-membered rings.
"""
from rdkit import Chem

def is_3_oxo_Delta_4__steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-Delta(4) steroid based on its SMILES string.
    The classification requires:
      (a) A fused steroid nucleus – defined by finding a connected set of rings (only considering rings
          of sizes 5 or 6) that comprises at least 4 rings and contains at least 15 carbons.
      (b) Within one of the six-membered rings of that fused core, an enone motif is present.
          To be precise, we require a substructure [C;R](=O)[C;R]=[C;R] found completely within the ring.
    Args:
        smiles (str): SMILES string of the molecule
    Returns:
        bool: True if the molecule is a 3-oxo-Delta(4) steroid, False otherwise.
        str: Explanation of the classification.
    """
    # Parse the SMILES into a molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # ---------------------------
    # 1. Identify candidate rings (only 5- or 6-membered rings are of interest)
    # ---------------------------
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()  # each is a tuple of atom indices
    candidate_rings = []
    for ring in all_rings:
        if len(ring) in (5, 6):
            candidate_rings.append(set(ring))
    if not candidate_rings:
        return False, "No rings of size 5 or 6 detected"

    # ---------------------------
    # 2. Build a fused ring graph.
    # We consider rings fused if they share at least one bond (i.e. at least 2 atoms).
    # ---------------------------
    num = len(candidate_rings)
    ring_graph = {i: set() for i in range(num)}
    for i in range(num):
        for j in range(i+1, num):
            # if rings share at least 2 atoms, consider them fused
            if len(candidate_rings[i] & candidate_rings[j]) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
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
                for nbr in ring_graph[current]:
                    if nbr not in comp:
                        stack.append(nbr)
            visited |= comp
            fused_components.append(comp)
    
    # ---------------------------
    # 3. Look for a candidate steroid nucleus.
    # We require:
    #    - At least 4 fused rings (the tetracyclic core),
    #    - Among which at least 3 are six-membered and at least 1 is five-membered,
    #    - And the union of atoms involved contains at least 15 carbon atoms.
    # ---------------------------
    steroid_component = None
    for comp in fused_components:
        comp_rings = [candidate_rings[i] for i in comp]
        if len(comp_rings) < 4:
            continue
        count6 = sum(1 for ring in comp_rings if len(ring)==6)
        count5 = sum(1 for ring in comp_rings if len(ring)==5)
        if count6 < 3 or count5 < 1:
            continue
        # build union of atoms in the connected component
        comp_atoms = set()
        for ring in comp_rings:
            comp_atoms.update(ring)
        # count carbons among these atoms
        c_count = sum(1 for idx in comp_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        if c_count < 15:
            continue
        steroid_component = comp_rings  # candidate nucleus obtained
        break

    if steroid_component is None:
        return False, "Steroid nucleus not found (insufficient fused ring system with 3 six-membered and 1 five-membered rings, or too few carbons)"
    
    # Build a union set of atom indices in the steroid nucleus.
    steroid_atoms = set()
    for ring in steroid_component:
        steroid_atoms.update(ring)
    
    # ---------------------------
    # 4. Look for an enone (α,β-unsaturated ketone) motif in one of the six-membered rings.
    # We want to identify a pattern: a carbon (in-ring) double-bonded to an oxygen,
    # which is connected (via a single bond) to another carbon that in turn is double-bonded to a third carbon.
    #
    # We define a SMARTS that represents C(=O)C=C all in a ring.
    # Note: [C;R] restricts the atom to be in a ring.
    # ---------------------------
    enone_smarts = "[C;R](=O)[C;R]=[C;R]"
    enone_query = Chem.MolFromSmarts(enone_smarts)
    
    # Only check six-membered rings from the steroid nucleus.
    six_membered_rings = [ring for ring in steroid_component if len(ring)==6]
    enone_found = False
    for ring in six_membered_rings:
        # For each match of the enone pattern in the full molecule,
        # check if ALL atoms in the match lie within the ring.
        for match in mol.GetSubstructMatches(enone_query):
            if set(match).issubset(ring):
                enone_found = True
                break
        if enone_found:
            break

    if not enone_found:
        return False, "Enone motif not found in any six-membered ring of the steroid nucleus"

    return True, "Contains a fused steroid nucleus (with at least 4 rings and sufficient carbon content) and a 3-oxo/Δ(4) enone motif within a six-membered ring."

# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with one of the examples: (20S)-20-hydroxypregn-4-en-3-one
    test_smiles = "[H][C@@]1(CC[C@@]2([H])[C@]3([H])CCC4=CC(=O)CC[C@]4(C)[C@@]3([H])CC[C@]12C)[C@H](C)O"
    result, reason = is_3_oxo_Delta_4__steroid(test_smiles)
    print(result, reason)