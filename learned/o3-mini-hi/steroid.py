"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroid compounds based on a modified interpretation of the cyclopenta[a]phenanthrene nucleus.
Definition: Any compound (natural or synthetic) that is based on a cyclopenta[a]phenanthrene backbone or its analogue.
The algorithm looks for a fused ring cluster composed of non‐aromatic rings (sizes 5,6,7) that either form:
 • an intact steroid nucleus: at least 4 fused rings including exactly one 5-membered ring and at least 2 rings of size 6 or 7, 
   and most of the atoms in the fused cluster are carbon.
 • a “seco” steroid nucleus: a fused cluster of 3 rings including at least one 5-membered ring and covering at least 16 unique atoms,
   with most atoms being carbon.
Also requires a minimal overall carbon count.
"""

from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    
    The method verifies:
      - The molecule contains at least 17 carbon atoms.
      - There is at least one fused cluster of non‐aromatic rings of size 5–7.
      - Among such clusters, either an intact steroid nucleus (4 or more fused rings with exactly one 5‐membered and
        at least two 6–7‐membered rings) is present OR a partially broken (seco) nucleus with 3 fused rings, containing at least one 5‐membered ring and spanning at least 16 unique atoms.
      - In addition, the fused cluster must be composed predominantly (≥80%) of carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is classified as a steroid, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall carbon count.
    c_total = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_total < 17:
        return False, f"Too few carbons ({c_total}) for a steroid nucleus"
    
    # Get ring information.
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No ring system found"
    
    # Filter rings: keep only rings with allowed sizes (5,6,7) and ignore rings that are completely aromatic.
    allowed_sizes = {5, 6, 7}
    valid_rings = []
    for ring in all_rings:
        ring_size = len(ring)
        if ring_size not in allowed_sizes:
            continue
        # Check if the ring is not fully aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        valid_rings.append(ring)
    if not valid_rings:
        return False, "No non-aromatic rings of size 5-7 found"
    
    # Build a graph where each node is a ring (by index in valid_rings) and an edge exists if rings share at least 2 atoms.
    ring_graph = {i: set() for i in range(len(valid_rings))}
    for i in range(len(valid_rings)):
        for j in range(i+1, len(valid_rings)):
            if len(set(valid_rings[i]).intersection(valid_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components (fused ring clusters) using depth-first search (DFS).
    def dfs(node, visited, comp):
        visited.add(node)
        comp.add(node)
        for neighbor in ring_graph[node]:
            if neighbor not in visited:
                dfs(neighbor, visited, comp)
    
    visited = set()
    fused_components = []
    for i in ring_graph:
        if i not in visited:
            comp = set()
            dfs(i, visited, comp)
            fused_components.append(comp)
    
    # Evaluate each fused-ring component.
    for comp in fused_components:
        comp_rings = [valid_rings[i] for i in comp]
        ring_sizes = [len(r) for r in comp_rings]
        num_rings = len(comp_rings)
        count_5 = sum(1 for size in ring_sizes if size == 5)
        count_6_7 = sum(1 for size in ring_sizes if size in {6,7})
        # Get the unique atoms (set of indices) in this fused cluster.
        unique_atoms = set()
        for ring in comp_rings:
            unique_atoms.update(ring)
        num_unique = len(unique_atoms)
        
        # Calculate fraction of carbons in the fused cluster.
        carbon_in_cluster = sum(1 for idx in unique_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6)
        frac_carbon = carbon_in_cluster / num_unique if num_unique > 0 else 0
        
        # Only accept clusters that are predominantly carbon (at least 80%).
        if frac_carbon < 0.80:
            continue
        
        # Criterion 1: Intact steroid nucleus (preferably 4 or more fused rings).
        if num_rings >= 4:
            # For an intact nucleus, we require exactly one 5-membered ring (typically the D ring)
            # and at least two rings of size 6 (or 7).
            if count_5 == 1 and (num_rings - count_5) >= 2:
                return True, ("Contains an intact fused steroid nucleus: "
                              f"{num_rings} fused rings (ring sizes: {ring_sizes}), with 1 five‐membered and {num_rings - 1} six/seven‐membered rings "
                              f"covering {num_unique} atoms (carbon fraction: {frac_carbon:.2f}).")
        # Criterion 2: Partially opened (seco) steroid nucleus.
        if num_rings == 3:
            if count_5 >= 1 and count_6_7 >= 1 and num_unique >= 16:
                return True, ("Contains a partially opened (seco) steroid nucleus: "
                              f"{num_rings} fused rings (ring sizes: {ring_sizes}) covering {num_unique} atoms "
                              f"(carbon fraction: {frac_carbon:.2f}).")
    
    return False, "No steroid nucleus pattern detected"

# For simple testing:
if __name__ == "__main__":
    test_examples = [
        # True positives
        ("O[C@]12[C@]([C@H](OC(=O)/C=C(/C(C)C)\\C)C[C@]3([C@@]1(O)CC=C4[C@@]3(CC[C@H](O)C4)C)[H])([C@](O)(CC2)C(=O)C)C", "Caudatin"),
        ("O[C@@H]1[C@]2(C(C3C([C@@]4(C(CC3O)=CC(=O)CC4)C)C1)CC[C@@]2([C@@H](CCCC(C)C(O)=O)C)[H])C", "7alpha,12alpha-dihydroxy-3-oxocholest-4-en-27-oic acid"),
        # Secosteroid example
        ("O[C@]1([C@]2([C@@]([C@](CC2)([C@@H](CCCC(O)(C)C)C)[H])(CCC1)C)[H])/C=C/C=3C(CCC(=O)C3)=C", "Secosteroid example"),
        # False-positive (non-steroid fused system)
        ("C[C@H](CCC(=O)NCCCC(O)=O)[C@H]1CC[C@H]2[C@H]3[C@H](CC[C@]12C)[C@@]4(C)CC[C@@H](O)C[C@H]4CC3=O", "7-oxo-gamma-aminoisobutyrolithocholic acid"),
    ]
    
    for smi, name in test_examples:
        res, reason = is_steroid(smi)
        print(f"Test {name}: {res} -- {reason}")