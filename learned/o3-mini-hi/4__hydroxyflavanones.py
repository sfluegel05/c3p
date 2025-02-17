"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 4'-hydroxyflavanones
Definition: “Any hydroxyflavanone having a hydroxy substituent located at position 4'.”
Heuristic (improved): 
  (1) The molecule must contain a cyclic ketone that is part of a ring system that also contains at least one oxygen atom.
      This is used as a proxy for the flavanone (benzopyranone) core.
  (2) There must be at least one six‐membered aromatic ring (candidate B ring) that is NOT fused with the core 
      but is attached by a single bond to any atom in the core.
  (3) In that candidate B ring, one of the ring atoms (the attachment point) must have an -OH substituent located
      at the para position (i.e. 3 bonds away along the ring) relative to the attachment.
      
Any molecule passing these tests is classified as a 4'-hydroxyflavanone.
"""
from rdkit import Chem

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
    atom_rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    
    # Step 1. Identify the flavanone core:
    # We require a cyclic ketone (C=O) where the carbonyl carbon is in a ring that also contains at least one oxygen.
    core_atoms = None
    ketone_smarts = Chem.MolFromSmarts("[#6;R]=O")  # looks for a carbonyl on a ring; note: the O is not in the ring
    matches = mol.GetSubstructMatches(ketone_smarts)
    for match in matches:
        # match is a tuple with the carbonyl carbon index only.
        keto_c = match[0]
        # Look through rings that contain the ketone carbon.
        for ring in atom_rings:
            if keto_c in ring:
                # Check if this ring contains an oxygen atom (other than the carbonyl oxygen which is not in the ring)
                ring_has_oxygen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
                if ring_has_oxygen:
                    core_atoms = set(ring)
                    break
        if core_atoms is not None:
            break
    if core_atoms is None:
        return False, "No flavanone core detected (cyclic ketone in an oxygen-containing ring not found)"
    
    # Step 2. Now look for a candidate B ring.
    # Determine all six-membered rings that are fully aromatic.
    candidate_B_rings = []
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # not six-membered
        # Check if every atom in ring is aromatic:
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Exclude rings fused with the core (we want the B ring to be connected by a single bond).
        if core_atoms & set(ring):  # if any atom is in common, then this is fused so skip
            continue
        # Also, verify that the ring is attached to the core by a single bond:
        attached = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look at neighbors that are not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms:
                    attached = True
                    break
            if attached:
                candidate_B_rings.append(tuple(ring))
                break
    if not candidate_B_rings:
        return False, "No six-membered aromatic B ring attached to the flavanone core detected"
    
    # Helper: Build an adjacency dictionary (graph) for a given ring (list of atom indices)
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
    
    # Helper: Compute the shortest path length between two nodes in a small graph using BFS.
    def ring_distance(g, start, end):
        from collections import deque
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
        return None  # should not happen in a connected cycle
    
    # Step 3. For each candidate B ring, check if there is an -OH group attached at the para position relative to the attachment.
    for ring in candidate_B_rings:
        # Build an adjacency graph for the ring.
        g = build_ring_graph(ring)
        # Identify attachment atoms: those in the candidate ring that have a neighbor in the core.
        attachment_atoms = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms:
                    attachment_atoms.append(idx)
                    break
        if not attachment_atoms:
            continue
        
        # For each atom in the ring, check if it bears an -OH group.
        # We record such atoms (within the ring) along with their indices in the ring. 
        # (Since the ring tuple order is arbitrary, we use the graph distances on the ring.)
        hydroxy_atoms = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # Look among neighbors (outside the ring) for an oxygen that likely is -OH.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms in the same ring
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    hydroxy_atoms.append(idx)
                    break
        
        if not hydroxy_atoms:
            continue
        
        # Now, for each attachment atom, check if any hydroxy atom in the ring is para to it.
        # In a six-membered ring the para position is 3 bonds away along the ring graph.
        for attach_idx in attachment_atoms:
            for oh_idx in hydroxy_atoms:
                dist = ring_distance(g, attach_idx, oh_idx)
                if dist == 3:
                    return True, ("Molecule has a flavanone core (cyclic ketone in an oxygen heterocycle) and a six-membered aromatic B ring "
                                  "attached via a single bond with an -OH substituent para to the attachment (position 4').")
    return False, "No appropriate B ring with a para -OH substituent detected on the flavanone core"

# (Optional) Testing examples.
if __name__ == '__main__':
    test_smiles = [
        "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",  # (+)-taxifolin, should be True
        "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",             # sakuranetin, should be True
        "Oc1ccc(cc1)C1CC(=O)c2c(O)cc(O)cc2O1",                  # naringenin, should be True
        "COC(=O)C1=CC=CC=C1"                                   # An ester, not a flavanone, should be False
    ]
    for sm in test_smiles:
        result, reason = is_4__hydroxyflavanones(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")