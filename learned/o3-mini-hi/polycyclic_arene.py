"""
Classifies: CHEBI:33848 polycyclic arene
"""
#!/usr/bin/env python
"""
Classifies: Polycyclic Arene, defined as a polycyclic aromatic hydrocarbon.
A polycyclic arene is a molecule that contains two or more candidate (i.e. carbon‐only) aromatic rings 
that are fused either directly (sharing a common bond) or indirectly (bridged via a non‐candidate ring).
This should (mostly) cover cases such as naphthalene, anthracene, fluorene, and larger polycyclic 
hydrocarbons while filtering out aromatic heterocycles.
"""

from rdkit import Chem

def is_polycyclic_arene(smiles: str):
    """
    Determines if a molecule is a polycyclic arene (a polycyclic aromatic hydrocarbon)
    based on its SMILES string.
    
    The program works by:
      1. Parsing the SMILES string and obtaining the ring information.
      2. Identifying candidate rings that are fully aromatic AND consist solely of carbon atoms.
      3. Building a graph of all rings (using a strict fusion criterion, i.e. sharing at least 2 atoms).
      4. For each connected component of rings, counting how many candidate rings occur.
         If two or more candidate rings occur in one connected component, the molecule is 
         classified as a polycyclic arene.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if molecule is classified as polycyclic arene,
                     otherwise False. The second element is a string explaining the classification.
    """
    # Parse SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = [set(ring) for ring in ring_info.AtomRings()]  # each ring represented as set of atom indices
    if not atom_rings:
        return False, "No rings detected in the molecule"
    
    # Define a helper to decide if a ring is a candidate aromatic hydrocarbon ring.
    # We require that every atom in the ring is aromatic and is carbon (atomic number 6).
    def is_aromatic_carbon_ring(ring):
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if not atom.GetIsAromatic() or atom.GetAtomicNum() != 6:
                return False
        return True
    
    # Mark candidate rings (fully aromatic, carbon-only)
    candidate_flags = [is_aromatic_carbon_ring(ring) for ring in atom_rings]
    
    # Build a connectivity graph of all rings.
    # Two rings are considered "fused" if they share at least 2 atoms (i.e. a common bond).
    n = len(atom_rings)
    ring_graph = {i: set() for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            if len(atom_rings[i].intersection(atom_rings[j])) >= 2:
                ring_graph[i].add(j)
                ring_graph[j].add(i)
    
    # Find connected components in the ring graph.
    visited = set()
    for i in range(n):
        if i in visited:
            continue
        # DFS to get all rings connected to ring i
        stack = [i]
        component = set()
        while stack:
            node = stack.pop()
            if node not in component:
                component.add(node)
                stack.extend(ring_graph[node] - component)
        visited |= component
        
        # Count candidate rings in this connected component
        candidate_count = sum(1 for idx in component if candidate_flags[idx])
        if candidate_count >= 2:
            return True, f"Found a fused aromatic hydrocarbon system with {candidate_count} candidate rings"
    
    return False, "No fused aromatic hydrocarbon system detected"

# For testing purposes you could add:
if __name__ == '__main__':
    # Example test cases (expected True for known polycyclic arenes)
    test_smiles = [
        "c1ccc2c(c1)ccc1cc3c(ccc4ccccc34)cc21",  # dibenz[a,h]anthracene
        "c1ccc2cc3ccccc3cc2c1",                  # anthracene
        "C1C2=CC=CC=C2C2=CC=CC=C12"               # fluorene (note: may be challenging)
    ]
    for smi in test_smiles:
        result, reason = is_polycyclic_arene(smi)
        print(smi, "=>", result, ":", reason)