"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
Definition: Any 3β-hydroxy-steroid that contains a double bond between positions 5 and 6.
Improvement: In addition to the previous heuristics (a six-membered ring with a double bond and 
a ring-bound chiral carbon carrying an –OH) we try to ensure that these features come from a fused 
ring system of at least three rings (i.e. a steroid nucleus).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if the given SMILES string represents a 3β-hydroxy-Δ5-steroid.

    The method checks three conditions:
      1. The molecule must contain a fused steroid nucleus – here we require a cluster of at least 
         three rings (fused meaning they share at least two atoms).
      2. Within the fused nucleus, there must be at least one six-membered ring that contains a 
         C=C double bond (a proxy for the Δ5 feature).
      3. Also within the fused nucleus, there must be at least one chiral carbon that is part of 
         a ring and has an attached hydroxyl group (a proxy for the 3β-hydroxy substituent).

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Brief reason for classification.
    """

    # Parse the molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # (A) Collect ring atom sets from the molecule.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings found in molecule"
    
    # (B) Build a graph where each node is a ring (represented by its atom set)
    # Two rings are considered "fused" if they share at least 2 atoms.
    n_rings = len(ring_info)
    adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(ring_info[i])
        for j in range(i+1, n_rings):
            set_j = set(ring_info[j])
            if len(set_i.intersection(set_j)) >= 2:
                adj[i].add(j)
                adj[j].add(i)
    
    # (C) Find connected clusters (i.e. fused ring systems) using DFS.
    visited = set()
    clusters = []
    for i in range(n_rings):
        if i not in visited:
            cluster = set()
            stack = [i]
            while stack:
                cur = stack.pop()
                if cur in cluster:
                    continue
                cluster.add(cur)
                for nb in adj[cur]:
                    if nb not in cluster:
                        stack.append(nb)
            clusters.append(cluster)
            visited |= cluster
    
    # We now choose the largest fused ring cluster. For a steroid nucleus we require at least three rings.
    largest_cluster = max(clusters, key=lambda c: len(c))
    if len(largest_cluster) < 3:
        return False, "No fused ring system (steroid nucleus) of at least 3 rings found"
    
    # Create a set of atom indices part of the fused steroid nucleus.
    nucleus_atoms = set()
    for ring_idx in largest_cluster:
        nucleus_atoms.update(ring_info[ring_idx])
    
    # (D) Heuristic 1: Within the nucleus, find at least one six-membered ring that contains a double bond.
    delta5_found = False
    for ring_idx in largest_cluster:
        ring = ring_info[ring_idx]
        if len(ring) == 6:
            # For each bond in the ring, check if it is a double bond.
            ring_atoms = list(ring)
            for i in range(len(ring_atoms)):
                a1 = ring_atoms[i]
                a2 = ring_atoms[(i+1)%len(ring_atoms)]
                bond = mol.GetBondBetweenAtoms(a1, a2)
                if bond is not None and bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Also require that both atoms are in the fused nucleus (should normally be true)
                    if a1 in nucleus_atoms and a2 in nucleus_atoms:
                        delta5_found = True
                        break
            if delta5_found:
                break
    if not delta5_found:
        return False, "No six-membered ring with a double bond (Δ5 feature) found in the steroid nucleus"
    
    # (E) Heuristic 2: Within the nucleus, look for a ring-bound chiral carbon that carries an –OH.
    hydroxy_found = False
    for atom in mol.GetAtoms():
        # Only consider atoms that are in the nucleus
        if atom.GetIdx() not in nucleus_atoms:
            continue
        if atom.GetAtomicNum() == 6 and atom.IsInRing():
            # Check if chirality has been assigned (this is a rough proxy for beta orientation)
            if atom.HasProp('_CIPCode'):
                # Check neighbors for an oxygen that is part of an –OH (at least one hydrogen)
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8:
                        if neighbor.GetTotalNumHs() >= 1:
                            hydroxy_found = True
                            break
        if hydroxy_found:
            break
    if not hydroxy_found:
        return False, "No ring-bound chiral carbon with an attached hydroxyl group (3β-hydroxy) found in the nucleus"
    
    return True, ("Molecule contains a fused steroid nucleus (≥3 fused rings), a six-membered ring with a "
                  "double bond (Δ5 feature) and a ring-bound chiral carbon with a hydroxyl (proxy for 3β-hydroxy), "
                  "consistent with a 3β-hydroxy-Δ5-steroid.")