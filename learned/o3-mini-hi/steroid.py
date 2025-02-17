"""
Classifies: CHEBI:35341 steroid
"""
"""
Classifies: Steroid compounds based on a cyclopenta[a]phenanthrene‐derived nucleus.
Definition: Any compound (natural or synthetic) based on the cyclopenta[a]phenanthrene skeleton,
which traditionally contains three fused six-membered rings and one five-membered ring, but allowing
for bond scissions, ring expansions, and contractions.
"""

from rdkit import Chem

def is_steroid(smiles: str):
    """
    Determines if a molecule is a steroid based on its SMILES string.
    
    First, it tries a relaxed SMARTS for an intact cyclopenta[a]phenanthrene nucleus.
    If not found, it builds a fused ring graph: two rings are considered fused if they
    share at least two atoms. Then, for each connected fused‐ring cluster, we count as “ideal”
    those rings whose sizes are in {5,6,7} (allowing for slight modifications) and that are not
    fully aromatic. A molecule is considered steroid‐like if it has either:
      - A fused cluster of 3 rings with at least 2 ideal rings, or
      - A fused cluster of 4 or more rings with at least 3 ideal rings.
    
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if molecule is a steroid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Try to match an intact but relaxed steroid nucleus SMARTS.
    # This pattern is a heuristic: we look for a fused tetracyclic system with a 5-membered ring.
    # (Note: The pattern is deliberately relaxed.)
    steroid_smarts = "[C;R1]1[C;R][C;R][C;R]2[C;R][C;R][C;R][C;R]3[C;R][C;R][C;R][C;R]4[C;R]([C;R]3)[C;R]2[C;R]1[C;R]4"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if steroid_pattern is not None and mol.HasSubstructMatch(steroid_pattern):
        return True, "Contains an intact cyclopenta[a]phenanthrene nucleus (relaxed SMARTS match)"
    
    # Next, analyze fused ring systems.
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # each ring is a tuple of atom indices
    if not rings:
        return False, "No ring system found"
    
    # Build a graph where each node is a ring (by index in 'rings'),
    # and an edge connects two rings if they share at least 2 atoms.
    ring_graph = {i: set() for i in range(len(rings))}
    for i in range(len(rings)):
        for j in range(i+1, len(rings)):
            if len(set(rings[i]).intersection(rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components in the ring graph.
    def dfs(node, visited, component):
        visited.add(node)
        component.add(node)
        for neighbor in ring_graph[node]:
            if neighbor not in visited:
                dfs(neighbor, visited, component)
    
    visited = set()
    components = []
    for node in ring_graph:
        if node not in visited:
            comp = set()
            dfs(node, visited, comp)
            components.append(comp)
    
    # Define allowed ring sizes (5-, 6-, and 7-membered rings) typical in steroid nuclei.
    allowed_sizes = {5, 6, 7}
    
    # Analyze each fused-ring component.
    for comp in components:
        comp_rings = [rings[i] for i in comp]
        ideal_count = 0  # count rings with size in allowed_sizes and not fully aromatic
        for ring in comp_rings:
            ring_size = len(ring)
            if ring_size in allowed_sizes:
                # Check if ring is fully aromatic.
                # If all atoms in the ring are aromatic, we do not count it as ideal.
                if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    ideal_count += 1
        
        comp_size = len(comp)
        # Heuristic rules:
        # - For a three-ring system, require at least two ideal (non-aromatic, allowed size) rings.
        # - For a four-or-more ring system, require at least three ideal rings.
        if comp_size == 3 and ideal_count >= 2:
            return True, ("Contains a fused ring system resembling a steroid nucleus "
                          "(three fused rings with at least two rings of allowed size and non-aromatic)")
        if comp_size >= 4 and ideal_count >= 3:
            return True, ("Contains a fused ring system resembling a steroid nucleus "
                          "(four or more fused rings with at least three rings of allowed size and non-aromatic)")
    
    return False, "No steroid nucleus pattern detected"
    
# Example usage (for testing purposes):
if __name__ == "__main__":
    # Test with some known SMILES examples:
    steroid_smiles = "OC12C(C(CC1)C(OC3OC(C(O)C(O)C3O)CO)C)(CCC4C2CC=C5C4(CCC(OC6OC(C(OC7OC(C(O)C(O)C7O)CO)C(OC)C6O)C)C5)C)C"  # Russelioside B
    nonsteroid_smiles = "COC(=O)\\C=C/NC(=O)c1cc2c3ccccc3[nH]c2c(n1)[C@H](C)OC(=O)C(\\C)=C/C"  # Dichotomide X (false positive previously)
    
    result, reason = is_steroid(steroid_smiles)
    print("Steroid test:", result, reason)
    
    result, reason = is_steroid(nonsteroid_smiles)
    print("Non-steroid test:", result, reason)