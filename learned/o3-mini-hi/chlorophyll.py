"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: Chlorophyll – a family of magnesium porphyrins defined by the presence 
of a fifth ring (an extra isocyclic ring) beyond the four pyrrole-like rings and usually 
including a long phytol (aliphatic) chain.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is a chlorophyll based on its SMILES string.
    
    Chlorophyll is defined (approximately) as a magnesium porphyrin that meets the 
    following criteria:
      1. The molecule is neutral.
      2. Contains at least one magnesium atom.
      3. At least one magnesium atom is coordinated to four nitrogen atoms 
         (as in the porphyrin core).
      4. Contains at least five rings in total, and at least one of these rings is five-membered.
      5. Contains a long aliphatic (phytol) chain, here defined as a contiguous chain 
         of at least 10 carbon atoms outside any ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as chlorophyll, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the overall molecule is neutral by summing all formal charges.
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, f"Molecule net charge is {total_charge}, expected neutral"
    
    # Check for the presence of at least one magnesium (Mg, atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"
    
    # Verify that at least one Mg is coordinated to at least four nitrogen atoms.
    mg_with_four_n = False
    for mg in mg_atoms:
        n_neighbors = sum(1 for nbr in mg.GetNeighbors() if nbr.GetAtomicNum() == 7)
        if n_neighbors >= 4:
            mg_with_four_n = True
            break
    if not mg_with_four_n:
        return False, "Magnesium is not coordinated to at least four nitrogen atoms (porphyrin core missing)"
    
    # Check ring system. We expect at least 5 rings with one being 5-membered.
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        return False, f"Insufficient number of rings found ({num_rings}); expected at least 5"
    
    atom_rings = ring_info.AtomRings()
    has_five_membered = any(len(ring) == 5 for ring in atom_rings)
    if not has_five_membered:
        return False, "No five-membered ring (extra isocyclic ring) found"
    
    # Determine the set of atoms that belong to any ring.
    ring_atoms = set()
    for ring in atom_rings:
        ring_atoms.update(ring)
    
    # Build the subgraph of non-ring carbon atoms.
    # We allow a chain to start at a ring junction but then continue only through non-ring carbons.
    nonring_carbons = set()
    for idx, atom in enumerate(mol.GetAtoms()):
        # We focus on carbon atoms.
        if atom.GetAtomicNum() == 6:
            nonring_carbons.add(idx)
    # Note: Some carbons in rings may be part of a substituent chain if attached to non-ring carbons.
    # For our chain search, we want to traverse only atoms that are definitely in acyclic chains.
    # So we restrict to those carbon atoms that are not in a ring.
    acyclic_carbons = {idx for idx in nonring_carbons if idx not in ring_atoms}
    
    # Build adjacency list for the subgraph consisting of acyclic carbon atoms.
    carbon_graph = {}
    for idx in acyclic_carbons:
        carbon_graph[idx] = set()
    
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        # consider bond if both atoms are in the acyclic carbon set.
        if i in acyclic_carbons and j in acyclic_carbons:
            carbon_graph[i].add(j)
            carbon_graph[j].add(i)
    
    # Compute the longest simple path in the subgraph.
    # Because molecules are typically small, a DFS search is acceptable.
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

    if longest_chain < 10:
        return False, f"Longest acyclic carbon chain is only {longest_chain} atoms; expected at least 10 (phytol chain not found)"
    
    return True, ("Contains magnesium porphyrin core (Mg coordinated to >=4 N), extra five-membered ring and "
                  f"a long phytol chain (longest acyclic carbon chain = {longest_chain} atoms)")
    
# Example usage (can be commented out in production)
if __name__ == "__main__":
    # Example: chlorophyll a (should return True)
    test_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"
    result, reason = is_chlorophyll(test_smiles)
    print(result, reason)