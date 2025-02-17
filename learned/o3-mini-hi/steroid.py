"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroid compounds based on a modified interpretation of the cyclopenta[a]phenanthrene nucleus.
Definition: Any compound (natural or synthetic) that is based on a cyclopenta[a]phenanthrene (or analogue) backbone.
The algorithm relies on identifying a fused ring cluster meeting criteria for an intact steroid nucleus (4 rings),
or a “seco” (opened) nucleus (3 fused rings), plus the additional check for angular methyl substituents.
Note: This is an approximate method.
"""

from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    
    Our method checks:
      1. The molecule contains at least 17 carbons.
      2. It contains non‐aromatic rings (sizes 5, 6 or 7) that form a fused cluster
         (i.e. sharing at least 2 atoms per edge).
      3. If a fused cluster has 4 rings (intact nucleus), exactly one must be 5‑membered and the others 6–7‑membered,
         with at least 17 unique atoms and predominately carbon (fraction ≥ 0.90).
      4. If a cluster has 3 rings (a possible “seco” steroid), it must span at least 16 unique atoms.
      5. Additionally, we check for one or more angular CH3 groups (a carbon atom attached directly to a cluster atom,
         that is itself terminal) as a further indicator of a steroid nucleus.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a steroid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that substituent counts are more reliable.
    mol = Chem.AddHs(mol)
    
    # Check overall carbon count.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 17:
        return False, f"Too few overall carbons ({total_carbons}) for a steroid nucleus"
    
    # Grab ring information from the molecule.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Filter for rings with allowed sizes (5,6,7) that are not fully aromatic.
    allowed_sizes = {5, 6, 7}
    valid_rings = []
    for ring in all_rings:
        ring_size = len(ring)
        if ring_size not in allowed_sizes:
            continue
        # Skip rings that are completely aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        valid_rings.append(ring)
    if not valid_rings:
        return False, "No non-aromatic rings of size 5-7 found"
    
    # Build a graph where each node is a ring. Two rings are connected if they share at least 2 atoms.
    ring_graph = {i: set() for i in range(len(valid_rings))}
    for i in range(len(valid_rings)):
        for j in range(i+1, len(valid_rings)):
            if len(set(valid_rings[i]).intersection(valid_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Use DFS to get the connected fused ring clusters.
    def dfs(node, visited, comp):
        visited.add(node)
        comp.add(node)
        for nbr in ring_graph[node]:
            if nbr not in visited:
                dfs(nbr, visited, comp)
    
    visited = set()
    fused_clusters = []
    for i in ring_graph:
        if i not in visited:
            comp = set()
            dfs(i, visited, comp)
            fused_clusters.append(comp)
    
    # Define a helper function to count angular methyl groups.
    def count_angular_methyl(cluster_atom_idxs):
        count = 0
        for idx in cluster_atom_idxs:
            atom = mol.GetAtomWithIdx(idx)
            # Check neighbors that are not part of the fused cluster.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in cluster_atom_idxs:
                    continue
                # Look for a terminal methyl: carbon, degree=1 (only connected to this atom)
                if nbr.GetAtomicNum() == 6 and nbr.GetDegree() == 1:
                    count += 1
        return count

    # Evaluate each fused-cluster candidate.
    for comp in fused_clusters:
        comp_rings = [valid_rings[i] for i in comp]
        ring_sizes = [len(r) for r in comp_rings]
        num_rings = len(comp_rings)
        count_5 = sum(1 for size in ring_sizes if size == 5)
        count_6_7 = sum(1 for size in ring_sizes if size in {6, 7})
        
        # Collect unique atoms in this cluster.
        unique_atoms = set()
        for ring in comp_rings:
            unique_atoms.update(ring)
        num_unique = len(unique_atoms)
        
        # Calculate fraction of atoms (in the fused cluster) that are carbon.
        carbons_in_cluster = sum(1 for idx in unique_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        frac_carbon = carbons_in_cluster / num_unique if num_unique > 0 else 0
        
        # Count angular methyl substituents (simplest CH3 groups off the cluster).
        methyl_count = count_angular_methyl(unique_atoms)
        
        # Check intact steroid nucleus pattern: 4 fused rings, exactly one 5-membered and at least two 6/7 rings,
        # at least 17 unique atoms, high carbon fraction and at least one angular methyl.
        if num_rings == 4 and count_5 == 1 and (num_rings - count_5) >= 2 and num_unique >= 17 and frac_carbon >= 0.90 and methyl_count >= 1:
            return True, (f"Contains an intact fused steroid nucleus: {num_rings} fused rings (ring sizes: {ring_sizes}), "
                          f"with 1 five‐membered and {num_rings - 1} six/seven‐membered rings covering {num_unique} atoms "
                          f"(carbon fraction: {frac_carbon:.2f}), and {methyl_count} angular methyl substituent(s).")
        
        # Check for a secosteroid (partially opened nucleus): 3 fused rings, at least one 5-membered,
        # spanning at least 16 atoms and adequate carbon enrichment.
        if num_rings == 3 and count_5 >= 1 and num_unique >= 16 and frac_carbon >= 0.85:
            return True, (f"Contains a partially opened (seco) steroid nucleus: {num_rings} fused rings (ring sizes: {ring_sizes}), "
                          f"covering {num_unique} atoms (carbon fraction: {frac_carbon:.2f}).")
    
    return False, "No steroid nucleus pattern detected"

# For simple testing:
if __name__ == "__main__":
    test_examples = [
        # True positives (examples taken from our outcomes)
        ("O[C@]12[C@]([C@H](OC(=O)/C=C(/C(C)C)C)C[C@]3([C@@]1(O)CC=C4[C@@]3(CC[C@H](O)C4)C)[H])([C@](O)(CC2)C(=O)C)", "Caudatin"),
        ("O[C@@H]1[C@]2(C(C3C([C@@]4(C(CC3O)=CC(=O)CC4)C)C1)CC[C@@]2([C@@H](CCCC(C)C(O)=O)C)[H])C", "7alpha,12alpha-dihydroxy-3-oxocholest-4-en-27-oic acid"),
        ("OC12C(C(CC1)C(OC3OC(C(O)C(O)C3O)CO)C)(CCC4C2CC=C5C4(CCC(OC6OC(C(OC7OC(C(O)C(O)C7O)CO)C(OC)C6O)C)C5)C)C", "Russelioside B"),
        # Example secosteroid candidate: note that our method might classify these based on the 3-ring cluster.
        ("O[C@]1([C@]2([C@@]([C@](CC2)([C@@H](CCCC(O)(C)C)C)[H])(CCC1)C)[H])/C=C/C=3C(CCC(=O)C3)=C", "Secosteroid example"),
        # False positive candidate – our stricter criteria should avoid mis‐classifying these.
        ("C[C@H](CCC(=O)NCCCC(O)=O)[C@H]1CC[C@H]2[C@H]3[C@H](CC[C@]12C)[C@@]4(C)CC[C@@H](O)C[C@H]4CC3=O", "7-oxo-gamma-aminoisobutyrolithocholic acid"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_steroid(smi)
        print(f"Test {name}: {result} -- {reason}")