"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
#!/usr/bin/env python
"""
Classifies: 4'-hydroxyflavanones
Definition: “Any hydroxyflavanone having a hydroxy substituent located at position 4'.”
Revised approach:
  (1) Identify a flavanone core by searching for a cyclic ketone (C=O) that belongs to a six‐membered ring
      which also contains an oxygen and sufficient aromatic character (at least 3 aromatic atoms) as expected for 
      a fused benzopyranone.
  (2) Identify six‐membered aromatic rings (candidate B rings) that are not fused into the core but attached via a single 
      single bond from an sp³ (non‐aromatic) core atom.
  (3) Within such a B ring, check that the atom (in the B ring) directly attached to the core has an –OH substituent 
      in its para position (i.e. three bonds away around the ring).
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
       bool: True if the molecule is classified as a 4'-hydroxyflavanone, False otherwise.
       str: A reason explaining the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get info on rings. Each ring is represented as a tuple of atom indices.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings detected"

    # STEP 1: Identify the flavanone core.
    # Heuristic: The flavanone core should contain a cyclic ketone (C=O) in a 6-membered ring 
    # that also includes at least one oxygen in the ring and has enough aromatic atoms (>=3).
    core_atoms = None
    ketone_smarts = Chem.MolFromSmarts("[#6;R]=O")  # a carbonyl where the carbon is in a ring
    matches = mol.GetSubstructMatches(ketone_smarts)
    for match in matches:
        keto_c = match[0]
        # Look within rings that contain the ketone carbon; require ring to have exactly 6 atoms.
        for ring in atom_rings:
            if keto_c in ring and len(ring)==6:
                # Check that the ring has at least one oxygen and enough aromatic atoms.
                has_oxygen = any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 for idx in ring)
                aromatic_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetIsAromatic())
                if not has_oxygen:
                    continue
                if aromatic_count < 3:
                    continue
                core_atoms = set(ring)
                break
        if core_atoms is not None:
            break
    if core_atoms is None:
        return False, ("No appropriate flavanone core was detected – the molecule lacks a 6-membered ring containing "
                       "a ketone, an internal oxygen and sufficient aromatic character (fused benzene ring) typical of flavanones.")

    # STEP 2: Identify candidate B rings.
    # Look for six-membered aromatic rings that do not share atoms with the core.
    candidate_B_rings = []
    for ring in atom_rings:
        if len(ring) != 6:
            continue
        if any(idx in core_atoms for idx in ring):
            continue  # skip rings fused with the core
        # Must be fully aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        
        # Count bonds connecting the ring to the core.
        connecting_bonds = 0
        core_attachment_atom = None
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                # if a neighbor is in the core then we count an inter-ring connection
                if nbr.GetIdx() in core_atoms:
                    connecting_bonds += 1
                    core_attachment_atom = nbr  # keep the last encountered core neighbor
        # require exactly one connection
        if connecting_bonds != 1:
            continue
        # In typical flavanones the core atom connected to the B ring should be sp3 (i.e. not aromatic)
        if core_attachment_atom is None or core_attachment_atom.GetIsAromatic():
            continue
        candidate_B_rings.append(tuple(ring))
    if not candidate_B_rings:
        return False, ("No appropriate six-membered aromatic B ring was detected that is attached via a single bond "
                       "from a non‐aromatic (sp³) atom in the core.")

    # Helper: Build an undirected graph (adjacency list) for atoms in a given ring.
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

    # Helper: Compute shortest path distance within the ring graph between two atoms (using BFS).
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
        return None  # should not happen in a connected ring

    # STEP 3: In each candidate B ring, verify that the connection point on the B ring
    # carries an –OH substituent at the para position (3 bonds away along the ring).
    for ring in candidate_B_rings:
        g = build_ring_graph(ring)
        # Identify the attachment atom in the B ring: it is the one that has a neighbor in the core.
        attachment_atoms = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in core_atoms:
                    attachment_atoms.append(idx)
                    break
        # There should be exactly one point of attachment.
        if len(attachment_atoms) != 1:
            continue
        attach_idx = attachment_atoms[0]
        
        # In the B ring, search for an atom that carries a hydroxyl (-OH) group 
        # (i.e. has an external oxygen with at least one hydrogen), 
        # and check that it is exactly 3 bonds away from the attachment atom (para position).
        found_para_OH = False
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # search through neighbors that are not part of the ring
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # require the substituent to be O and (by typical RDKit count) to have at least one hydrogen
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    # Check ring distance from attachment atom to the current atom (within the ring)
                    d = ring_distance(g, attach_idx, idx)
                    if d == 3:
                        found_para_OH = True
                        break
            if found_para_OH:
                break
        if found_para_OH:
            return True, ("Molecule has a flavanone core (a 6‐membered ring containing a ketone, an oxygen and sufficient aromatic character) "
                          "and a non‐fused six‐membered aromatic B ring attached via a single bond from an sp³ core atom "
                          "with a hydroxyl substituent at the para (4') position.")
    
    return False, "No candidate B ring with a para –OH substituent (at 4' position) was detected on the flavanone core."


# (Optional) Testing examples.
if __name__ == '__main__':
    test_smiles = [
        "O[C@@H]1[C@H](Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",  # (+)-taxifolin - expected True
        "COc1cc(O)c2C(=O)C[C@H](Oc2c1)c1ccc(O)cc1",            # (2S)-2'-methoxykurarinone - True
        "Oc1ccc2C(=O)C[C@H](Oc2c1)c1ccc(O)c(O)c1",            # butin - True
        "O1[C@@H](CC(=O)C=2C1=C(CC=C(C)C)C(O)=CC2O)C3=CC(=C(O)C=C3)CC=C(C)C",  # Euchrestaflavanone A - false positive in previous attempt
        "COC(=O)C1=CC=CC=C1"                                  # not a flavanone - False
    ]
    for sm in test_smiles:
        result, reason = is_4__hydroxyflavanones(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")