"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: Chlorophyll – a family of magnesium porphyrins (or chlorins)
defined by the presence of a fused macrocycle with at least five rings (the extra ring being the “fifth ring” beyond the four pyrrole‐like rings).
Although many chlorophylls feature a long phytol chain (or exocyclic ester/carboxyl group), some forms (chlorophyllides) lack it.
The heuristic therefore checks for (a) a magnesium atom, (b) a fused macrocycle with at least 5 rings (including one 5‐membered ring),
(c) exactly 4 nitrogen atoms in rings, and (d) reports side‐chain features if present.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is classified as chlorophyll (or a chlorophyllide) based on its SMILES string.
    
    The heuristic uses these steps:
      1. The molecule is parsed.
      2. The molecule must contain at least one magnesium atom.
      3. The molecule must feature a fused ring system of at least 5 rings and at least one ring of size 5.
      4. Within these rings, exactly 4 nitrogen atoms are found (the tetrapyrrole core).
      5. Check the magnesium bonding environment: if any Mg atom is bonded to 4 nitrogen atoms, record that;
         if not, we still accept but note that the Mg coordination is not ideal.
      6. Report side‐chain features: compute the longest connected acyclic (non‐ring) carbon chain and 
         search for exocyclic ester ([OX2;!R]-C(=O)[#6]) or carboxylic acid ([CX3](=O)[OX2H]) groups.
         These are not strictly required so as to accept chlorophyllides (which may lack a long phytol chain).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as chlorophyll/chlorophyllide, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for magnesium (Mg, atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"
    
    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # tuple of tuples of atom indices in each ring
    num_rings = len(atom_rings)
    if num_rings < 5:
        return False, f"Insufficient rings ({num_rings}); expected at least 5 in the fused chlorin core"
    
    # Check that at least one ring in the set is five-membered
    has_five_membered = any(len(ring) == 5 for ring in atom_rings)
    if not has_five_membered:
        return False, "No five-membered ring found in the fused system (expected as the extra ring)"
    
    # Count nitrogen atoms that are in any ring – expect exactly 4 (the tetrapyrrole core)
    ring_nitrogen_idxs = set()
    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                ring_nitrogen_idxs.add(idx)
    if len(ring_nitrogen_idxs) != 4:
        return False, f"Found {len(ring_nitrogen_idxs)} nitrogen atoms in rings; expected exactly 4 for the tetrapyrrole core"
    
    # Check magnesium coordination: look for any Mg atom directly bonded to nitrogen(s)
    mg_coordination = None
    for mg in mg_atoms:
        n_neighbors = sum(1 for nbr in mg.GetNeighbors() if nbr.GetAtomicNum() == 7)
        if n_neighbors == 4:
            mg_coordination = "Mg directly coordinated to 4 nitrogen atoms"
            break
    if mg_coordination is None:
        mg_coordination = "Mg present but not clearly coordinated to 4 nitrogen atoms"
    
    # Now report side-chain features.
    # Build set of indices of atoms that belong to any ring.
    ring_atom_idxs = set()
    for ring in atom_rings:
        ring_atom_idxs.update(ring)
    
    # Identify acyclic carbons (atomic number 6) not in any ring.
    acyclic_carbons = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in ring_atom_idxs}
    
    # Build a graph of connectivity among the acyclic carbons.
    carbon_graph = {idx: set() for idx in acyclic_carbons}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in acyclic_carbons and j in acyclic_carbons:
            carbon_graph[i].add(j)
            carbon_graph[j].add(i)
    
    # DFS search to determine the length of the longest acyclic carbon chain.
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
    
    # Search for exocyclic ester and free acid groups.
    ester_smarts = "[OX2;!R]-C(=O)[#6]"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_query)
    exo_ester_found = False
    for match in ester_matches:
        # Check if the oxygen is exocyclic: it is not in a ring, but at least one neighbor is in a ring.
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        if not o_atom.IsInRing() and any(nbr.IsInRing() for nbr in o_atom.GetNeighbors()):
            exo_ester_found = True
            break

    acid_smarts = "[CX3](=O)[OX2H]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    acid_found = False
    for match in acid_matches:
        c_idx = match[0]
        c_atom = mol.GetAtomWithIdx(c_idx)
        if not c_atom.IsInRing() and any(nbr.IsInRing() for nbr in c_atom.GetNeighbors()):
            acid_found = True
            break

    # Build a side-chain commentary but do not require a minimum chain length.
    side_chain_comments = []
    if longest_chain > 0:
        side_chain_comments.append(f"longest acyclic carbon chain is {longest_chain} atoms")
    if exo_ester_found:
        side_chain_comments.append("exocyclic ester group present")
    if acid_found:
        side_chain_comments.append("free carboxylic acid group present")
    side_chain_comment = "; ".join(side_chain_comments) if side_chain_comments else "no significant side-chain features found"
    
    # Summarize the findings.
    reason = ("Contains magnesium; " + mg_coordination +
              f"; fused macrocycle with {num_rings} rings (including a five-membered ring); " +
              f"found {len(ring_nitrogen_idxs)} ring nitrogens; side-chain features: {side_chain_comment}.")
    
    return True, reason

# For testing
if __name__ == "__main__":
    test_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"
    result, explanation = is_chlorophyll(test_smiles)
    print("Result:", result)
    print("Explanation:", explanation)