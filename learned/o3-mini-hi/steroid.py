"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroid compounds 
Definition: Any compound (natural or synthetic) whose core can be traced to a cyclopenta[a]phenanthrene backbone (or a related variant). 
This algorithm attempts to identify a fused ring cluster (3 or 4 rings) having a number of atoms in the range expected for a steroid nucleus, predominantly carbons, and showing one or more angular methyl substituents.
Note: This is an approximate method.
"""

from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.

    Our method:
      1. Parses the SMILES and adds explicit hydrogens.
      2. Looks for rings of sizes 5, 6 or 7 that are non-aromatic.
      3. Builds a "fused" ring graph from rings that share at least 2 atoms.
      4. Searches for fused ring clusters with at least 3 rings.
      5. For each fused cluster, we require:
           • The cluster must contain at least 16 unique atoms.
           • Most of the atoms (≥80%) are carbon.
           • At least one angular methyl (a terminal CH3 group) is attached to a cluster atom.
      6. If any cluster meets these criteria, we call the molecule a steroid.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a steroid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for more reliable substituent counts.
    mol = Chem.AddHs(mol)
    
    # Get overall ring information.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings found in the molecule"
    
    # Filter for rings with allowed sizes (5, 6 and 7) that are NOT fully aromatic.
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
        return False, "No non‐aromatic rings of size 5–7 found"

    # Build a graph where each node represents a ring.
    # Two rings are connected if they share at least 2 atom indices.
    ring_graph = {i: set() for i in range(len(valid_rings))}
    for i in range(len(valid_rings)):
        for j in range(i+1, len(valid_rings)):
            if len(set(valid_rings[i]).intersection(valid_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Use depth-first search (DFS) to detect connected fused ring clusters.
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
    
    # Helper function to count angular (terminal) methyl groups attached to atoms in the cluster.
    def count_angular_methyl(cluster_atom_idxs):
        count = 0
        for idx in cluster_atom_idxs:
            atom = mol.GetAtomWithIdx(idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in cluster_atom_idxs:
                    continue
                # Accept terminal (degree 1) carbon substituents.
                if nbr.GetAtomicNum() == 6 and nbr.GetDegree() == 1:
                    count += 1
        return count

    # Evaluate each fused ring cluster.
    for comp in fused_clusters:
        # comp is a set of indices of valid_rings that are fused.
        comp_rings = [valid_rings[i] for i in comp]
        # Gather unique atoms that are part of the fused rings.
        unique_atoms = set()
        for ring in comp_rings:
            unique_atoms.update(ring)
        num_unique = len(unique_atoms)
        
        # Count how many atoms in the cluster are carbon.
        carbon_count = sum(1 for idx in unique_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        frac_carbon = carbon_count / num_unique if num_unique else 0
        
        # Count the number of rings in this cluster.
        num_rings = len(comp_rings)
        
        # Count angular methyl substituents.
        methyl_count = count_angular_methyl(unique_atoms)
        
        # For a steroid nucleus (or its opened "seco" variant), we expect:
        # - a fused ring cluster of at least 3 rings,
        # - covering roughly 16 or more atoms,
        # - a high carbon fraction (>= 0.80),
        # - and at least one terminal methyl group.
        if num_rings >= 3 and num_unique >= 16 and frac_carbon >= 0.80 and methyl_count >= 1:
            return True, (f"Contains a fused steroid‐like nucleus with {num_rings} rings covering {num_unique} atoms "
                          f"(carbon fraction: {frac_carbon:.2f}), with {methyl_count} angular methyl substituent(s).")
    
    return False, "No steroid nucleus pattern detected"

# For simple testing:
if __name__ == "__main__":
    test_examples = [
        # Example steroids: a few examples (names are as provided in the description):
        ("O=C(N[C@@H](C(O)(C)C)C)[C@H]([C@@]1([C@@]2([C@@](CC1)(/C(/CCC2)=C/C=C\\3/C[C@@H](O)C[C@H](O)C3=C)[H])C)[H])C", 
         "1alpha,25-dihydroxy-24-oxo-23-azavitamin D2 / 1alpha,25-dihydroxy-24-oxo-23-azaergocalciferol"),
        ("O[C@]12[C@]([C@H](OC(=O)/C=C(/C(C)C)C)C[C@]3([C@@]1(O)CC=C4[C@@]3(CC[C@H](O)C4)C)[H])([C@](O)(CC2)C(=O)C)", 
         "Caudatin"),
        ("O[C@@H]1[C@]2(C(C3C([C@@]4(C(CC3O)=CC(=O)CC4)C)C1)CC[C@@]2([C@@H](CCCC(C)C(O)=O)C)[H])C", 
         "7alpha,12alpha-dihydroxy-3-oxocholest-4-en-27-oic acid"),
        ("OC12C(C(CC1)C(OC3OC(C(O)C(O)C3O)CO)C)(CCC4C2CC=C5C4(CCC(OC6OC(C(OC7OC(C(O)C(O)C7O)CO)C(OC)C6O)C)C5)C)C", 
         "Russelioside B"),
    ]
    
    for smi, name in test_examples:
        result, reason = is_steroid(smi)
        print(f"Test {name}: {result} -- {reason}")