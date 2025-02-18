"""
Classifies: CHEBI:28966 chlorophyll
"""
"""
Classifies: Chlorophyll – a family of magnesium porphyrins (or chlorins) defined by the presence
of a fused macrocycle that contains a fifth (extra, non‐pyrrole) ring. The rings have various side chains,
most often including a long phytol chain (or exocyclic ester groups in chlorophyllides).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_chlorophyll(smiles: str):
    """
    Determines if a molecule is classified as chlorophyll (or chlorophyllide) based on its SMILES string.
    
    The decision uses the following heuristic:
      1. The molecule must be valid and neutral.
      2. The molecule must contain at least one magnesium atom.
      3. The overall fused ring system must comprise at least 5 rings, and at least one should be five-membered.
      4. Count the number of nitrogen atoms present in rings (i.e. part of the macrocycle).
         We require at least 4 such N atoms even if the “direct coordination” to Mg is not obvious.
      5. Finally, to filter out false positives (e.g. magnesium protoporphyrins), we check that the molecule has
         either a “long” acyclic carbon chain (≥10 C atoms outside rings, as a phytol chain) or it shows at least one
         exocyclic ester group ([OX2;!R]-C(=O)) attached to the macrocycle.
         
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as chlorophyll (or chlorophyllide), False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall neutrality (sum of formal charges should be zero)
    total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    if total_charge != 0:
        return False, f"Molecule net charge is {total_charge}, expected neutral"

    # Check for magnesium (atomic number 12)
    mg_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 12]
    if not mg_atoms:
        return False, "No magnesium atom found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 5:
        return False, f"Insufficient number of rings found ({num_rings}); expected at least 5"

    # Get list of all rings (as tuples of atom indices)
    atom_rings = ring_info.AtomRings()
    
    # Check that at least one ring is five-membered.
    if not any(len(ring) == 5 for ring in atom_rings):
        return False, "No five-membered ring (extra fused ring) found"

    # Gather nitrogen atoms that are part of any ring
    ring_nitrogens = set()
    for ring in atom_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 7:
                ring_nitrogens.add(idx)
    
    if len(ring_nitrogens) < 4:
        return False, f"Only {len(ring_nitrogens)} nitrogen atoms found in rings; expected at least 4"
    
    # OPTIONAL: now check for the presence of a characteristic side chain.
    # We define two mutually exclusive criteria:
    # (a) A long acyclic carbon chain (phytol chain) outside of any ring – we define long as 10 or more connected carbons.
    # (b) At least one exocyclic ester group attached (SMARTS: [OX2;!R]-C(=O)[#6]).
    # (A molecule meeting either (a) or (b) is acceptable.)
    
    # First, identify atoms that belong to any ring.
    ring_atom_idxs = set()
    for ring in atom_rings:
        ring_atom_idxs.update(ring)
    
    # Build an acyclic carbon graph of carbons not in rings
    acyclic_carbons = {atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIdx() not in ring_atom_idxs}
    carbon_graph = {idx: set() for idx in acyclic_carbons}
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        if i in acyclic_carbons and j in acyclic_carbons:
            carbon_graph[i].add(j)
            carbon_graph[j].add(i)
            
    # DFS-based search to find longest path in the acyclic carbon subgraph.
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

    # Check for exocyclic ester groups (attached to rings but not in rings).
    ester_smarts = "[OX2;!R]-C(=O)[#6]"
    ester_query = Chem.MolFromSmarts(ester_smarts)
    ester_matches = mol.GetSubstructMatches(ester_query)
    
    # Flag whether one of the ester atoms is attached to a ring atom (assume excocyclic if the O is not in a ring but neighbors a ring atom)
    exo_ester_found = False
    for match in ester_matches:
        # match is a tuple of atom indices for the pattern [O;!R]-C(=O)-...
        o_idx = match[0]
        o_atom = mol.GetAtomWithIdx(o_idx)
        # if any neighbor of this O is in a ring, we consider this ester exocyclic to a cyclic system.
        if any(nbr.GetIdx() in ring_atom_idxs for nbr in o_atom.GetNeighbors()):
            exo_ester_found = True
            break

    # Now decide: require at least one of the two criteria is met.
    chain_ok = longest_chain >= 10
    ester_ok = exo_ester_found

    if not (chain_ok or ester_ok):
        chain_comment = f"No long phytol chain detected (longest acyclic carbon chain is {longest_chain} atoms) and no exocyclic ester group found"
    else:
        if chain_ok:
            chain_comment = f"Contains long acyclic chain of {longest_chain} carbons (phytol chain present)"
        else:
            chain_comment = "No long acyclic chain detected but exocyclic ester group(s) present"

    # Construct a reason message that summarizes key features.
    reason = ("Contains magnesium; fused ring system with at least 5 rings (including a five-membered ring); "
              f"found {len(ring_nitrogens)} ring nitrogens; {chain_comment}.")
    
    return True, reason

# Example (this block can be commented out when used as a module)
if __name__ == "__main__":
    # Test a chlorophyll a structure (should be True)
    test_smiles = "CCC1=C(C)C2=Cc3c(C=C)c(C)c4C=C5[C@@H](C)[C@H](CCC(=O)OC\\C=C(/C)CCC[C@H](C)CCC[C@H](C)CCCC(C)C)C6=[N+]5[Mg--]5(n34)n3c(=CC1=[N+]25)c(C)c1C(=O)[C@H](C(=O)OC)C6=c31"
    result, reason = is_chlorophyll(test_smiles)
    print(result)
    print(reason)