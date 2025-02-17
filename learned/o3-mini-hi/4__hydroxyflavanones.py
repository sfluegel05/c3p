"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 4'-hydroxyflavanones
Definition: “Any hydroxyflavanone having a hydroxy substituent located at position 4'.”
Revised approach:
  (1) Identify a flavanone core by searching for a cyclic ketone (C=O) in a ring that also contains at least one oxygen,
      and importantly, the ring must be at least partially saturated (i.e. not all atoms in it are aromatic).
  (2) Identify six‐membered aromatic rings (candidate B rings) not fused with the core and attached by exactly one 
      single bond from an sp³ (non‐aromatic) atom in the core.
  (3) In the candidate B ring, check that the attachment point’s para position (atoms 3 bonds away along the ring) 
      bears an –OH group.
Any molecule passing all tests is classified as a 4'-hydroxyflavanone.
"""

from rdkit import Chem
from collections import deque

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if the molecule represented by `smiles` is a 4'-hydroxyflavanone.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if molecule is classified as a 4'-hydroxyflavanone, False otherwise.
       str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of atom index tuples

    # Step 1: Identify flavanone core.
    # Heuristic: Look for a cyclic ketone ([#6;R]=O) whose ring (the carbonyl-bearing ring) also has an oxygen
    # and is partly saturated (i.e. not all atoms aromatic).
    core_atoms = None
    ketone_smarts = Chem.MolFromSmarts("[#6;R]=O")
    matches = mol.GetSubstructMatches(ketone_smarts)
    for match in matches:
        keto_c = match[0]  # carbon in C=O
        # Look for a ring containing the ketone carbon.
        for ring in atom_rings:
            if keto_c in ring:
                # Check that the ring has at least one oxygen atom (could belong to the ring)
                has_oxygen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
                if not has_oxygen:
                    continue
                # Ensure that the ring is not entirely aromatic.
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    continue
                core_atoms = set(ring)
                break
        if core_atoms is not None:
            break
    if core_atoms is None:
        return False, "No appropriate flavanone core detected (cyclic ketone in a partly unsaturated, oxygen-containing ring not found)"
    
    # Step 2: Identify candidate B rings:
    # We look for six-membered rings that are fully aromatic and do not share atoms with the core.
    candidate_B_rings = []
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        if any(idx in core_atoms for idx in ring):
            continue  # fused with core; skip
        # Must be fully aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Count the number of bonds connecting this ring to the core.
        connecting_bonds = 0
        core_attachment_atom = None
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms:
                    connecting_bonds += 1
                    core_attachment_atom = nbr
        # We require exactly one connection.
        if connecting_bonds != 1:
            continue
        # Also, ensure that the neighbor in the core is sp3 (i.e. not aromatic) to favor the flavanone (dihydro) connectivity.
        if core_attachment_atom is not None and core_attachment_atom.GetIsAromatic():
            continue
        candidate_B_rings.append(tuple(ring))
    if not candidate_B_rings:
        return False, "No appropriate six-membered aromatic B ring attached via a single bond to the (non‐aromatic) core detected"
    
    # Helper: Build a graph (adjacency list) for the atoms in a ring.
    def build_ring_graph(ring):
        g = {i: set() for i in ring}
        ring_set = set(ring)
        for bond in mol.GetBonds():
            a1 = bond.GetBeginAtomIdx()
            a2 = bond.GetEndAtomIdx()
            if a1 in ring_set and a2 in ring_set:
                g[a1].add(a2)
                g[a2].add(a1)
        return g

    # Helper: Compute shortest path in the ring graph between two atoms using BFS.
    def ring_distance(g, start, end):
        visited = {start}
        queue = deque([(start, 0)])
        while queue:
            current, d = queue.popleft()
            if current == end:
                return d
            for nbr in g[current]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, d+1))
        return None  # Should not happen for a connected ring

    # Step 3: For each candidate B ring, verify that the atom of attachment on the ring
    # has a hydroxyl (-OH) substituent on the opposite (para) position (i.e. 3 bonds away in the ring).
    for ring in candidate_B_rings:
        g = build_ring_graph(ring)
        # Identify the attachment atom in the B ring: the one having a neighbor in the core.
        attachment_atoms = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms:
                    attachment_atoms.append(idx)
                    break
        if len(attachment_atoms) != 1:
            continue  # must be exactly one point of attachment
        attach_idx = attachment_atoms[0]
        
        # Identify atoms within the ring that bear a hydroxyl (-OH) substituent.
        hydroxy_positions = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look for external O with at least one H (typical for -OH)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # ignore atoms in the ring
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    hydroxy_positions.append(idx)
                    break
        if not hydroxy_positions:
            continue

        # In a six-membered ring, the para position is 3 bonds away.
        for oh_idx in hydroxy_positions:
            dist = ring_distance(g, attach_idx, oh_idx)
            if dist == 3:
                return True, ("Molecule has a flavanone core (cyclic ketone in a partly saturated, oxygen-containing ring) and a six-membered aromatic B ring "
                              "attached via a single bond (from a sp³ core atom) with an -OH substituent para to the attachment (position 4').")
    
    return False, "No candidate B ring with a para -OH substituent detected on the flavanone core"

# (Optional) Testing examples.
if __name__ == '__main__':
    test_smiles = [
        "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",  # (+)-taxifolin, expected: True
        "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",            # (2S)-2'-methoxykurarinone, True
        "Oc1ccc2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1",            # butin, True
        "COC(=O)C1=CC=CC=C1"                                  # an ester, not a flavanone, False
    ]
    for sm in test_smiles:
        result, reason = is_4__hydroxyflavanones(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")