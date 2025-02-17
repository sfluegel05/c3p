"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: CHEBI:3beta-hydroxy-Delta(5)-steroid
Definition: Any 3β-hydroxy-steroid that contains a double bond between positions 5 and 6.
Improvement Strategy:
  1. Parse the molecule (with stereochemistry and kekulization).
  2. Identify fused ring clusters (rings sharing ≥2 atoms) and select the largest cluster (≥3 rings) as the steroid nucleus.
  3. Collect six-membered rings from the nucleus.
     a. Mark rings that have exactly one nonaromatic double bond (candidate for the Δ5 bond location).
     b. Also, among all six-membered rings (including ones not bearing the double bond), search for a candidate ring that contains
        a chiral sp³ carbon bearing a “free” –OH (an oxygen with ≥1 hydrogen) and not participating in a double bond.
  4. In order for the molecule to be classified as a 3β-hydroxy-Δ5-steroid, we require that the ring harboring the free –OH is fused
     (shares at least one atom) with one of the rings that contains exactly one nonaromatic double bond.
  5. If such fused rings exist, the molecule is classified as a 3β-hydroxy-Δ5-steroid.
  
If any step cannot be met then the function returns False with the appropriate explanation.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Determines if a molecule (given its SMILES string) is a 3β-hydroxy-Δ5-steroid.
    
    The improved algorithm:
      1. Parses the SMILES, assigns stereochemistry and kekulizes.
      2. Finds fused ring clusters (rings sharing at least 2 atoms) and selects the largest cluster (with ≥3 rings) as the candidate steroid nucleus.
      3. Among the six-membered rings in the nucleus, two types are identified:
           - Double Bond Rings: a six-membered ring with exactly one nonaromatic double bond.
           - OH Rings: a six-membered ring that has at least one carbon atom that is:
                 • atomic number 6,
                 • sp³-hybridized,
                 • stereochemically defined (chiral tag not unspecified),
                 • not directly involved in any double bond,
             and which carries at least one attached oxygen atom that explicitly has one hydrogen (a free –OH).
      4. The steroid hypothesis is met if any OH ring is fused (shares at least one atom) with a double bond ring.
    
    Args:
      smiles (str): SMILES string representing the molecule.
    
    Returns:
      (bool, str): A tuple with True and a reason if classified as a 3β-hydroxy-Δ5-steroid,
                   otherwise False with an explanation.
    """
    # Step 1: Parse the molecule and assign stereochemistry.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass  # Proceed even if kekulization fails

    # Step 2: Get ring information as tuples of atom indices.
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings detected"
    
    # Build an adjacency between rings: two rings are fused if they share at least 2 atoms.
    n_rings = len(ring_info)
    ring_adj = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        set_i = set(ring_info[i])
        for j in range(i+1, n_rings):
            set_j = set(ring_info[j])
            if len(set_i.intersection(set_j)) >= 2:
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Use DFS to find fused clusters.
    visited = set()
    clusters = []
    for i in range(n_rings):
        if i in visited: 
            continue
        cluster = set()
        stack = [i]
        while stack:
            cur = stack.pop()
            if cur in cluster:
                continue
            cluster.add(cur)
            stack.extend(ring_adj[cur] - cluster)
        clusters.append(cluster)
        visited |= cluster

    # Choose the largest fused ring cluster as the candidate steroid nucleus.
    largest_cluster = max(clusters, key=lambda x: len(x))
    if len(largest_cluster) < 3:
        return False, "No fused ring system (steroid nucleus) of at least 3 rings found"

    # Collect all atom indices in the steroid nucleus.
    nucleus_atoms = set()
    for ring_idx in largest_cluster:
        nucleus_atoms.update(ring_info[ring_idx])
    
    # Helper: Build a ring graph for a given set of atom indices (only within that ring).
    def build_ring_graph(ring_atoms):
        graph = {a: set() for a in ring_atoms}
        for a in ring_atoms:
            atom = mol.GetAtomWithIdx(a)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring_atoms:
                    bond = mol.GetBondBetweenAtoms(a, nbr_idx)
                    if bond is not None:
                        graph[a].add(nbr_idx)
        return graph

    # Helper: Compute shortest path length in a subgraph (BFS).
    def shortest_path_length(graph, start, end):
        if start == end:
            return 0
        visited_nodes = {start}
        queue = [(start, 0)]
        while queue:
            current, dist = queue.pop(0)
            for neighbor in graph[current]:
                if neighbor == end:
                    return dist + 1
                if neighbor not in visited_nodes:
                    visited_nodes.add(neighbor)
                    queue.append((neighbor, dist+1))
        return float('inf')
    
    # Step 3: Among rings in the nucleus, focus on six-membered rings.
    six_membered = {}
    for idx in largest_cluster:
        ring = list(ring_info[idx])
        if len(ring) == 6:
            six_membered[idx] = ring
    
    if not six_membered:
        return False, "No six-membered rings found in the steroid nucleus"
    
    # Identify "double bond rings": six-membered rings that contain exactly one nonaromatic double bond.
    double_bond_rings = {}
    for ring_idx, ring in six_membered.items():
        # Build ring graph for this ring.
        ring_graph = build_ring_graph(ring)
        seen_pairs = set()
        dbl_count = 0
        dbl_bond_atoms = set()
        for a in ring:
            for b in ring_graph[a]:
                pair = tuple(sorted((a, b)))
                if pair in seen_pairs:
                    continue
                seen_pairs.add(pair)
                bond = mol.GetBondBetweenAtoms(a, b)
                if bond is None:
                    continue
                if bond.GetBondType() == Chem.BondType.DOUBLE and (not bond.GetIsAromatic()):
                    dbl_bond_atoms.update([a, b])
                    dbl_count += 1
        if dbl_count == 1:
            double_bond_rings[ring_idx] = {"atoms": ring, "dbl_atoms": dbl_bond_atoms}
    
    if not double_bond_rings:
        return False, "No six-membered ring in the nucleus with exactly one nonaromatic double bond (Δ5 candidate) found"
    
    # Identify "OH rings": six-membered rings that have at least one candidate chiral sp³ carbon (not in any double bond) with an attached free OH.
    oh_rings = {}
    for ring_idx, ring in six_membered.items():
        candidates = []
        for atom_idx in ring:
            # Skip if atom is directly part of any double bond in any ring (we will check later if needed)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                continue
            # Check that this atom is not directly in any double bond (in its bonds)
            in_double = False
            for bond in atom.GetBonds():
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    in_double = True
                    break
            if in_double:
                continue
            # Look for an attached oxygen with at least one hydrogen (free –OH)
            has_free_OH = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check explicit hydrogen count on the oxygen.
                    if nbr.GetTotalNumHs() >= 1:
                        has_free_OH = True
                        break
            if has_free_OH:
                candidates.append(atom_idx)
        if candidates:
            oh_rings[ring_idx] = candidates
    
    if not oh_rings:
        return False, "No six-membered ring in the nucleus with a candidate chiral sp³ carbon bearing a free –OH found"
    
    # Step 4: Require that one of the OH rings and one of the double bond rings are fused.
    # Fusion here means that the two rings share at least one atom.
    for dbl_ring_idx, dbl_data in double_bond_rings.items():
        dbl_atoms = set(dbl_data["atoms"])
        for oh_ring_idx, oh_candidates in oh_rings.items():
            # If the rings are not the same, check fusion.
            if dbl_ring_idx == oh_ring_idx:
                # In many steroids the OH (in position 3) is on a different ring than the Δ5 double bond (position 5-6);
                # if they occur in the same ring, skip this candidate.
                continue
            oh_atoms = set(six_membered[oh_ring_idx])
            if dbl_atoms.intersection(oh_atoms):
                # Found a fused pair: one ring with the double bond and a different ring with the free –OH candidate.
                return True, ("Molecule contains a fused steroid nucleus (≥3 fused rings) with a six-membered ring "
                              "bearing exactly one nonaromatic double bond and an adjacent six-membered ring harboring a chiral "
                              "sp³ carbon with a free –OH group. This pattern is consistent with a 3β-hydroxy-Δ5-steroid.")
    
    return False, ("No appropriate fusion found between a six-membered ring with a nonaromatic double bond and a six-membered ring "
                   "carrying a chiral sp³ carbon with a free –OH group. This structural motif is required for a 3β-hydroxy-Δ5-steroid.")

# Example usage:
if __name__ == "__main__":
    # Example SMILES strings taken from the provided list.
    test_smiles = [
        # Expected positive (for instance, cholesterol-like):
        "C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",
        # A representative false positive candidate from the previous attempt:
        "CC[C@H](C)C(=O)[C@@H]1[C@H]2[C@@H]3[C@@H](O)OC(C)=CC3=CC(=O)[C@]2(C)OC1=O"
    ]
    for s in test_smiles:
        flag, reason = is_3beta_hydroxy_Delta_5__steroid(s)
        print("SMILES:", s)
        print("Classification:", flag)
        print("Reason:", reason)
        print("------------------------------------------------")