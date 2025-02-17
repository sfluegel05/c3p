"""
Classifies: CHEBI:1722 3beta-hydroxy-Delta(5)-steroid
"""
"""
Classifies: 3beta-hydroxy-Delta(5)-steroid
Definition: Any 3β-hydroxy-steroid that contains a double bond between positions 5 and 6.
Improvement: 
  - The molecule is parsed with stereochemistry assignment and kekulized.
  - We first require a fused ring system of at least 3 rings (candidate steroid nucleus).
  - For each six-membered ring in that nucleus we require exactly one nonaromatic double bond.
  - Then we search that ring for a candidate sp³ carbon that is chiral and carries a “free” –OH (an attached oxygen having at least one hydrogen).
  - Instead of using the ring order we build a graph (only among atoms in the ring) and compute the shortest path
    distance between the candidate carbon and the two atoms in the double bond. The candidate is only accepted if 
    both distances are >= 2.
If a ring is found with these features then the molecule is classified as a 3β-hydroxy-Δ5-steroid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3beta_hydroxy_Delta_5__steroid(smiles: str):
    """
    Checks whether the provided SMILES corresponds to a 3β-hydroxy-Δ5-steroid.
    
    The algorithm performs the following:
    1. Parses the SMILES, assigns stereochemistry and kekulizes the molecule.
    2. Extracts the ring information and groups rings into fused clusters (sharing at least 2 atoms).
       The largest fused cluster (with ≥3 rings) is treated as the candidate steroid nucleus.
    3. For each six-membered ring in that cluster:
       - The bonds within the ring are examined. We require exactly one nonaromatic double bond.
       - A graph is built among the ring atoms so that the “ring distance” (in bonds) can be computed.
       - For each atom (candidate carbon) in the ring that is:
            (a) a carbon (atomic number 6),
            (b) sp³-hybridized,
            (c) marked with a specific chiral tag (i.e. not CHI_UNSPECIFIED),
            (d) not directly part of the double bond,
         we check for an attached oxygen that has at least one hydrogen (a free –OH group).
         Then, using the ring-graph distance, we require that the candidate atom is at least two bonds
         away (along the ring) from each atom in the double bond.
    4. If at least one six-membered ring meets these conditions, the molecule is classified as a 3β-hydroxy-Δ5-steroid.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      bool: True if the molecule is classified as a 3β-hydroxy-Δ5-steroid, else False.
      str: Explanation for the classification decision.
    """
    # Step 1: Parse and preprocess the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass  # continue even if kekulization fails

    # Step 2: Get ring information (each ring is a tuple of atom indices)
    ring_info = mol.GetRingInfo().AtomRings()
    if not ring_info:
        return False, "No rings detected"

    # Build connectivity between rings: two rings are fused if they share at least 2 atoms.
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

    # Choose the largest fused ring cluster as the steroid nucleus candidate.
    largest_cluster = max(clusters, key=lambda x: len(x))
    if len(largest_cluster) < 3:
        return False, "No fused ring system (steroid nucleus) of at least 3 rings found"
    
    # Collect all atom indices in the steroid nucleus.
    nucleus_atoms = set()
    for ring_idx in largest_cluster:
        nucleus_atoms.update(ring_info[ring_idx])
    
    # Helper: For a given list of atom indices that form a ring, build a connectivity graph (dict) using bonds that are in the ring.
    def build_ring_graph(ring_atoms):
        graph = {a: set() for a in ring_atoms}
        for a in ring_atoms:
            atom = mol.GetAtomWithIdx(a)
            for nbr in atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx in ring_atoms:
                    # Check that there is a bond between these atoms
                    bond = mol.GetBondBetweenAtoms(a, nbr_idx)
                    if bond is not None:
                        graph[a].add(nbr_idx)
        return graph

    # Helper: Compute shortest path length between two atoms in a subgraph (BFS).
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
    
    # Step 3: For each six-membered ring in the nucleus candidate(s), try to find the characteristic motif.
    for ring_idx in largest_cluster:
        ring = list(ring_info[ring_idx])
        if len(ring) != 6:
            continue  # only consider six-membered rings
        
        # Build connectivity (a graph) for the atoms in this ring.
        ring_graph = build_ring_graph(ring)
        
        # Identify bonds within the ring that are nonaromatic double bonds.
        double_bond_atoms = set()
        double_bond_count = 0
        # Iterate over each bond between adjacent atoms (as found from ring_graph).
        # To avoid duplicates we check each pair only once.
        seen_pairs = set()
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
                    double_bond_atoms.update([a, b])
                    double_bond_count += 1
        # For a typical Δ5-steroid ring we expect exactly one nonaromatic double bond.
        if double_bond_count != 1:
            continue

        # For each atom in the ring, look for candidate: a chiral sp3 carbon not directly in the double bond.
        for atom_idx in ring:
            if atom_idx in double_bond_atoms:
                continue  # cannot be part of the double bond
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:
                continue
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                continue
            # Require explicit stereochemistry
            if atom.GetChiralTag() == Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                continue
            # Look for an attached oxygen that is “free” (has at least one hydrogen).
            has_free_OH = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # Check if oxygen bears an explicit hydrogen.
                    if nbr.GetTotalNumHs() >= 1:
                        has_free_OH = True
                        break
            if not has_free_OH:
                continue
            
            # Now, check the distance within the ring between this candidate and every atom involved in the double bond.
            meets_distance = True
            for db_atom in double_bond_atoms:
                d = shortest_path_length(ring_graph, atom_idx, db_atom)
                # Require at least 2 bonds along the ring
                if d < 2:
                    meets_distance = False
                    break
            if not meets_distance:
                continue
            
            # We found a six-membered ring with exactly one nonaromatic double bond and a candidate chiral sp3 carbon
            # with a free –OH located at least 2 bonds away (along the ring) from the double bond.
            return True, ("Molecule contains a fused steroid nucleus (≥3 fused rings) and at least one six-membered ring "
                          "with exactly one nonaromatic double bond and a chiral sp³ carbon bearing a free –OH group located "
                          "at least two bonds away (along the ring) from that double bond. This is consistent with a "
                          "3β-hydroxy-Δ5-steroid.")
    
    return False, ("No six-membered ring in the steroid nucleus was found that simultaneously has exactly one nonaromatic "
                   "double bond and a chiral sp³ carbon with a free –OH group (located at least two bonds away on the ring) "
                   "required for a 3β-hydroxy-Δ5-steroid.")
    
# Example usage:
if __name__ == "__main__":
    # Note: These examples come from the literature.
    test_smiles = [
        "C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C",  # (24S,25S)-cholest-5-en-3beta,24,26-triol (expected positive)
        "CC[C@H](C)C(=O)[C@@H]1[C@H]2[C@@H]3[C@@H](O)OC(C)=CC3=CC(=O)[C@]2(C)OC1=O"  # Example false positive candidate from previous failed attempt.
    ]
    for s in test_smiles:
        flag, reason = is_3beta_hydroxy_Delta_5__steroid(s)
        print("SMILES:", s)
        print("Classification:", flag)
        print("Reason:", reason)
        print("------------------------------------------------")