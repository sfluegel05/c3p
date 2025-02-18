"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: Chlorophyll – a family of magnesium porphyrins defined by the presence
of a fifth (non-pyrrole) ring beyond the four pyrrole-like rings. Although these molecules
usually include a long aliphatic (phytol) chain, its absence (as in chlorophyllides) is allowed.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is classified as chlorophyll based on its SMILES string.
    
    Chlorophyll (and chlorophyllides) are defined (approximately) as follows:
      1. The molecule is neutral.
      2. Contains at least one magnesium atom.
      3. At least one magnesium atom is coordinated to at least four nitrogen atoms (the porphyrin core).
      4. Has a fused ring system containing at least five rings overall, with at least one being five-membered 
         (the extra isocyclic ring beyond the four pyrrole-like rings).
      5. Ideally contains a long aliphatic chain (phytol chain), defined as a contiguous chain of at least 10 carbon atoms
         outside any ring. However, the absence of such a chain (as in chlorophyllides) will not cause rejection.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as chlorophyll (or chlorophyllide), False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check neutrality (sum of formal charges should be zero)
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, f"Molecule net charge is {total_charge}, expected neutral"
    
    # Check for the presence of magnesium (Mg, atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"
    
    # Verify that at least one Mg is coordinated to (roughly) 4 nitrogen atoms.
    # (Some chlorophyll-like compounds may show slight deviations; we require at least 4.)
    mg_with_four_n = False
    for mg in mg_atoms:
        n_neighbors = sum(1 for nbr in mg.GetNeighbors() if nbr.GetAtomicNum() == 7)
        if n_neighbors >= 4:
            mg_with_four_n = True
            break
    if not mg_with_four_n:
        return False, "Magnesium is not coordinated to at least four nitrogen atoms (porphyrin core missing)"
    
    # Check that the fused ring system has at least five rings in total.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        return False, f"Insufficient number of rings found ({num_rings}); expected at least 5"
    
    # Check for at least one five-membered ring among the ring systems
    atom_rings = ring_info.AtomRings()
    if not any(len(ring)==5 for ring in atom_rings):
        return False, "No five-membered ring (extra isocyclic ring) found"
    
    # Now, check for the presence of a long acyclic carbon chain (potential phytol chain).
    # Note: Chlorophyllides may lack the phytol chain so we use this as an optional flag.
    # Determine atoms that are in rings.
    ring_atoms = set()
    for ring in atom_rings:
        ring_atoms.update(ring)
    
    # Identify all carbon atoms
    carbon_indices = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6}
    # Consider only carbons not in any ring for chain search.
    acyclic_carbons = {idx for idx in carbon_indices if idx not in ring_atoms}
    
    # Build an adjacency graph for acyclic carbon atoms (by index).
    carbon_graph = {idx: set() for idx in acyclic_carbons}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in acyclic_carbons and j in acyclic_carbons:
            carbon_graph[i].add(j)
            carbon_graph[j].add(i)
    
    # Compute the longest simple path in this acyclic carbon subgraph using DFS.
    longest_chain = 0
    def dfs(current, visited):
        max_length = 1
        for neighbor in carbon_graph.get(current, []):
            if neighbor not in visited:
                length = 1 + dfs(neighbor, visited | {neighbor})
                if length > max_length:
                    max_length = length
        return max_length

    for start in carbon_graph:
        chain_length = dfs(start, {start})
        if chain_length > longest_chain:
            longest_chain = chain_length
    
    # The long aliphatic chain is optional. If present, we expect it to be at least 10 carbons.
    if longest_chain >= 10:
        chain_comment = f"Contains long acyclic chain of {longest_chain} carbons (phytol chain present)"
    else:
        chain_comment = f"No long phytol chain detected (longest acyclic carbon chain is {longest_chain} atoms)"
    
    # Return True with a reason that includes the key features.
    reason = ("Contains magnesium porphyrin core (Mg coordinated to >=4 N) and a fused ring system with at least 5 rings "
              "including a five-membered ring. " + chain_comment)
    return True, reason

# Example usage (this section can be commented out if not needed in production)
if __name__ == "__main__":
    # Test with chlorophyll a SMILES (should return True)
    test_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"
    result, reason = is_chlorophyll(test_smiles)
    print(result, reason)