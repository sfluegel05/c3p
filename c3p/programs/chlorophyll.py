"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: Chlorophyll – a family of magnesium porphyrins (or chlorins)
defined by a fused macrocycle (with at least five rings, one of which is five‐membered),
exactly four nitrogen atoms in rings, and usually a long phytol chain or a characteristic
exocyclic ester/carboxylic acid side group.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is classified as chlorophyll (or a chlorophyllide) based on its SMILES string.
    
    The heuristic uses these steps:
      1. The molecule must parse and be neutral.
      2. It must contain at least one magnesium atom.
      3. At least one magnesium atom must be directly connected to exactly 4 nitrogen atoms.
      4. The molecule must have a fused ring system of at least 5 rings, and at least one ring should be five-membered.
      5. Within the rings, exactly 4 nitrogen atoms must be present (the tetrapyrrole core).
      6. To avoid false positives (e.g. magnesium porphyrins or chlorophyllins), we require that
         the molecule shows one of the following side‐chain features:
             a. A long acyclic carbon chain (≥ 10 carbons) outside any ring,
             b. An exocyclic ester group (SMARTS: [OX2;!R]-C(=O)[#6]),
             c. Or a free carboxylic acid group (SMARTS: [CX3](=O)[OX2H]).
         (Many chlorophylls carry a long phytol chain; many chlorophyllides lack it but have a free acid.)
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if classified as chlorophyll/chlorophyllide, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall neutrality by summing formal charges
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, f"Molecule net charge is {total_charge}, expected neutral"
    
    # Check for magnesium (Mg atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"
    
    # Check if at least one Mg is directly connected to exactly 4 nitrogen atoms.
    mg_ok = False
    for mg in mg_atoms:
        # Count neighbors with atomic number 7 (nitrogen)
        n_neighbors = sum(1 for nbr in mg.GetNeighbors() if nbr.GetAtomicNum() == 7)
        if n_neighbors == 4:
            mg_ok = True
            break
    if not mg_ok:
        return False, "No magnesium atom directly coordinated to 4 nitrogen atoms"
    
    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        return False, f"Insufficient number of rings ({num_rings}); expected at least 5 in the fused system"
    
    atom_rings = ring_info.AtomRings()
    
    # Ensure at least one ring is five-membered
    if not any(len(ring) == 5 for ring in atom_rings):
        return False, "No five-membered ring (expected as extra fused ring) found"
    
    # Count nitrogen atoms that belong to rings (should equal 4 for a tetrapyrrole core)
    ring_nitrogens = set()
    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                ring_nitrogens.add(idx)
    if len(ring_nitrogens) != 4:
        return False, f"Found {len(ring_nitrogens)} nitrogen atoms in rings; expected exactly 4 for the tetrapyrrole core"
    
    # Now look for chlorophyll-like side chains.
    # First, build a set of indices for atoms in any ring.
    ring_atom_idxs = set()
    for ring in atom_rings:
        ring_atom_idxs.update(ring)
    
    # Identify acyclic carbon atoms (atomic number 6) not in any ring.
    acyclic_carbons = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in ring_atom_idxs}
    # Build an adjacency graph (only for these carbons)
    carbon_graph = {idx: set() for idx in acyclic_carbons}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in acyclic_carbons and j in acyclic_carbons:
            carbon_graph[i].add(j)
            carbon_graph[j].add(i)
    
    # DFS-based search to find the longest chain in the acyclic carbon subgraph.
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

    # Check for exocyclic ester groups.
    ester_smarts = "[OX2;!R]-C(=O)[#6]"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_query)
    exo_ester_found = False
    for match in ester_matches:
        # match[0] is the oxygen in the ester; if it is not in a ring but attached to a ring atom then count it.
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        if any(nbr.GetIdx() in ring_atom_idxs for nbr in o_atom.GetNeighbors()):
            exo_ester_found = True
            break

    # Check for free carboxylic acid group.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_query = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_query)
    acid_found = False
    # Only consider acid groups where the carbon is attached to a ring atom.
    for match in acid_matches:
        c_idx = match[0]
        c_atom = mol.GetAtomWithIdx(c_idx)
        if any(nbr.GetIdx() in ring_atom_idxs for nbr in c_atom.GetNeighbors()):
            acid_found = True
            break

    # Decide if side chain is acceptable.
    side_chain_ok = (longest_chain >= 10) or exo_ester_found or acid_found
    
    chain_comment = ""
    if longest_chain >= 10:
        chain_comment = f"Contains a long acyclic chain of {longest_chain} carbons (phytol chain present)"
    elif exo_ester_found:
        chain_comment = "Contains an exocyclic ester group attached to the macrocycle"
    elif acid_found:
        chain_comment = "Contains a free carboxylic acid group indicative of chlorophyllide"
    
    if not side_chain_ok:
        return False, (f"Side-chain filter not passed: longest acyclic carbon chain is {longest_chain} atoms, "
                       "and no exocyclic ester or free acid group found")
    
    # Summarize the findings.
    reason = ("Contains magnesium; at least one Mg is coordinated to 4 Ns; fused macrocycle with at least 5 rings "
              f"(including a five-membered ring); found 4 ring nitrogens; {chain_comment}.")
    
    return True, reason

# Example usage:
if __name__ == "__main__":
    # You can test a known chlorophyll/chlorophyllide SMILES here.
    test_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"
    result, explanation = is_chlorophyll(test_smiles)
    print("Result:", result)
    print("Explanation:", explanation)