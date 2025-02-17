"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroid compounds based on a modified interpretation of the cyclopenta[a]phenanthrene nucleus.
Definition: Any compound (natural or synthetic) that is based on the cyclopenta[a]phenanthrene backbone,
or a derived analogue of it (including secosteroids), is expected to have a fused ring system of non‐aromatic rings of size 5–7.
We require that at a minimum the molecule has at least 17 carbons and one fused‐ring cluster that either contains:
  - Four rings (an intact nucleus) with one 5‐membered ring and at least two 6–7‐membered rings; or
  - Three fused rings (a partially broken or “seco” nucleus) that include one 5‐membered ring and cover at least 16 unique atoms.
"""

from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    
    The method first checks that the molecule has the minimum number of carbon atoms
    expected for a steroid (>16). Then it extracts all rings from the molecule and
    only considers rings that are not fully aromatic and that have 5, 6, or 7 atoms (allowing for slight modifications).
    A ring–graph is constructed where two rings are considered fused if they share at least 2 atoms.
    Connected components (fused-ring clusters) are examined and if one component meets either
    of two criteria it is considered steroid-like:
      - An intact steroid nucleus: a cluster with 4 fused rings that contains exactly one 5-membered ring 
        (the typical D ring) and at least two rings of size 6 or 7.
      - A partially opened (seco) steroid: a cluster with 3 fused rings containing one 5-membered ring and 
        covering at least 16 unique atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as steroid, False otherwise.
        str: Explanation of the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall carbon count: steroids generally have many carbons.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 17:
        return False, f"Too few carbons ({c_count}) for a steroid nucleus"

    # Get ring information – each ring is a tuple of atom indices.
    ring_info = mol.GetRingInfo()
    complete_rings = ring_info.AtomRings()
    if not complete_rings:
        return False, "No ring system found"
    
    # Filter rings: only count those rings which do not consist solely of aromatic atoms
    # and whose sizes are typical for a steroid nucleus (5,6,7).
    allowed_sizes = {5, 6, 7}
    valid_rings = []
    for ring in complete_rings:
        ring_size = len(ring)
        if ring_size not in allowed_sizes:
            continue
        # Only accept the ring if not every atom is aromatic.
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        valid_rings.append(ring)
    if not valid_rings:
        return False, "No non-aromatic rings of size 5-7 found"
    
    # Build a graph connecting rings that share at least 2 atoms (fused rings).
    ring_graph = {i: set() for i in range(len(valid_rings))}
    for i in range(len(valid_rings)):
        for j in range(i+1, len(valid_rings)):
            if len(set(valid_rings[i]).intersection(valid_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components within the ring graph.
    def dfs(node, visited, comp):
        visited.add(node)
        comp.add(node)
        for neighbor in ring_graph[node]:
            if neighbor not in visited:
                dfs(neighbor, visited, comp)
    
    visited = set()
    components = []
    for i in ring_graph:
        if i not in visited:
            comp = set()
            dfs(i, visited, comp)
            components.append(comp)
    
    # Examine each fused-ring component.
    for comp in components:
        comp_rings = [valid_rings[i] for i in comp]
        ring_sizes = [len(ring) for ring in comp_rings]
        # Count how many rings are 5-membered.
        count_5 = sum(1 for size in ring_sizes if size == 5)
        count_6_7 = sum(1 for size in ring_sizes if size in {6, 7})
        comp_size = len(comp_rings)
        
        # Also count unique atoms in the component (for partially broken steroid nuclei)
        unique_atoms = set()
        for ring in comp_rings:
            unique_atoms.update(ring)
        
        # Criterion 1: intact steroid nucleus (four fused rings: one 5-membered and at least two 6/7 rings).
        if comp_size == 4:
            if count_5 == 1 and count_6_7 >= 2:
                return True, ("Contains an intact fused steroid nucleus: 4 non-aromatic rings "
                              f"(ring sizes: {ring_sizes}), including a single 5-membered ring.")
        # Criterion 2: a partially opened steroid nucleus (secosteroid) with 3 fused rings 
        # that include a 5-membered ring and cover enough atoms.
        if comp_size == 3:
            if count_5 >= 1 and count_6_7 >= 1 and len(unique_atoms) >= 16:
                return True, ("Contains a partially opened steroid nucleus: 3 fused non-aromatic rings "
                              f"(ring sizes: {ring_sizes}) covering {len(unique_atoms)} atoms.")
    
    return False, "No steroid nucleus pattern detected"


# For simple testing:
if __name__ == "__main__":
    # Example test cases (from provided outcomes):
    examples = [
        # A known steroid-like molecule (Caudatin)
        ("O[C@]12[C@]([C@H](OC(=O)/C=C(/C(C)C)\\C)C[C@]3([C@@]1(O)CC=C4[C@@]3(CC[C@H](O)C4)C)[H])([C@](O)(CC2)C(=O)C)C", "Caudatin"),
        # A known secosteroid: (6E)-(8S)-8,25-dihydroxy-9,10-seco-4,6,10(19)-cholestatrien-3-one (should be steroid even though ring opened)
        ("O[C@]1([C@]2([C@@]([C@](CC2)([C@@H](CCCC(O)(C)C)C)[H])(CCC1)C)[H])/C=C/C=3C(CCC(=O)C3)=C", 
         "Secosteroid example"),
        # A false positive example from previous attempt.
        ("C1[C@@H](O[C@@H]([C@H]2[C@@H]1C3=C(O2)C=CC(=C3)NS(=O)(=O)C4=CC=CC=C4)CO)CC(=O)NCC(F)(F)F", 
         "Non-steroid false positive"),
    ]
    
    for smi, name in examples:
        result, reason = is_steroid(smi)
        print(f"Test {name}: {result} -- {reason}")