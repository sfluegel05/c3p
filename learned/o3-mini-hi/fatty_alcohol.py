"""
Classifies: CHEBI:24026 fatty alcohol
"""
#!/usr/bin/env python
"""
Classifies: fatty alcohol
Definition: "An aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
             Fatty alcohols may be saturated or unsaturated and may be branched or unbranched."
             
This implementation:
  1. Parses the molecule using RDKit.
  2. Searches for candidate alcohol groups on non‐ring carbons with a directly attached hydroxyl (-OH).
  3. For each candidate alcohol carbon, it builds a connectivity graph of non‐aromatic, non‐ring carbons.
  4. Instead of using DFS to pick one branch, a BFS is used from the alcohol carbon to determine the maximum number
     of contiguous carbons reachable along any path.
  5. If any candidate alcohol group gives a longest chain (distance) of at least MIN_CHAIN_LENGTH carbons, the molecule
     is classified as a fatty alcohol.
     
Note:
  This method allows branches and both saturated and unsaturated chains.
"""

from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule belongs to the fatty alcohol class.

    Args:
      smiles (str): SMILES string of the molecule.

    Returns:
      (bool, str): A tuple; True if classified as a fatty alcohol, False otherwise,
                   plus a human‐readable reason.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS query for an alcohol group on a non‐ring carbon.
    # This captures a carbon (not in a ring) connected to a hydroxyl group.
    alcohol_smarts = "[C;!R][OX2H]"
    alcohol_query = Chem.MolFromSmarts(alcohol_smarts)
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    if not alcohol_matches:
        return False, "No proper alcohol (-OH) group found on a non‐ring carbon"
    
    # Build a set of carbon atoms that are non‐aromatic and not in rings.
    carbon_nodes = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and (not atom.IsInRing()):
            carbon_nodes.add(atom.GetIdx())
            
    # Build a connectivity graph among these carbon atoms (using indices).
    carbon_graph = {idx: [] for idx in carbon_nodes}
    for idx in carbon_nodes:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in carbon_nodes:
                carbon_graph[idx].append(nbr.GetIdx())
    
    # For each candidate alcohol carbon, use BFS to determine the furthest distance (number of atoms)
    # reachable along the connected acyclic aliphatic carbon network.
    # (This gives us the maximum contiguous chain length starting from the candidate alcohol carbon.)
    def bfs_max_chain_length(start):
        from collections import deque
        visited = {start}
        queue = deque([(start, 1)])  # (node, distance) distance counts the starting carbon as 1.
        max_length = 1
        while queue:
            current, dist = queue.popleft()
            if dist > max_length:
                max_length = dist
            # Visit neighbors from the carbon graph.
            for nbr in carbon_graph.get(current, []):
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, dist + 1))
        return max_length
    
    # Minimum chain length requirement.
    # Definition says "a chain of 3 to greater than 27 carbon atoms" so use 3 as the minimum.
    MIN_CHAIN_LENGTH = 3
    
    candidate_found = False
    best_chain_length = 0
    candidate_reason = ""
    
    # Iterate over candidate alcohol groups.
    for match in alcohol_matches:
        alc_c_idx = match[0]   # candidate alcohol carbon
        alc_o_idx = match[1]   # hydroxyl oxygen
        
        alc_atom = mol.GetAtomWithIdx(alc_c_idx)
        # Skip candidate if the alcohol carbon is in a ring (should be aliphatic).
        if alc_atom.IsInRing():
            continue
        
        # Check that apart from the -OH, all other neighbors of the alcohol carbon are carbons.
        valid = True
        for nbr in alc_atom.GetNeighbors():
            if nbr.GetIdx() == alc_o_idx:
                continue  # this is the hydroxyl oxygen
            if nbr.GetAtomicNum() != 6:
                valid = False
                break
        if not valid:
            continue
        
        # Ensure the candidate alcohol carbon is in our aliphatic carbon graph.
        if alc_c_idx not in carbon_graph:
            continue
        
        # Compute the maximum chain length from the candidate alcohol carbon using BFS.
        chain_length = bfs_max_chain_length(alc_c_idx)
        if chain_length > best_chain_length:
            best_chain_length = chain_length
        
        # We do not exclude branches since branching is allowed for fatty alcohols.
        if chain_length >= MIN_CHAIN_LENGTH:
            candidate_found = True
            candidate_reason = (f"Molecule contains an alcohol group with a contiguous acyclic aliphatic chain "
                                f"of {chain_length} carbons (starting at atom index {alc_c_idx}).")
            break  # Accept the first candidate meeting the criteria.
    
    if candidate_found:
        return True, candidate_reason
    else:
        return False, (f"Longest acyclic aliphatic chain from any candidate alcohol group is only "
                       f"{best_chain_length} carbons (need at least {MIN_CHAIN_LENGTH}).")

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCC(C)O",       # octan-2-ol
        "O[C@H](CCCCCCCCCCC(C)C)C",  # (S)-13-Methyltetradecan-2-ol
        "CC(C)CCCCCCCCCCCCCCCCCO",   # 18-methylicosan-1-ol
        "OCCCCCCCCCCCCCCCC#CCC=C",    # 17-OCTADECEN-14-YN-1-OL
        "CCCCCCCCCC(O)CCCCCCCCCCC",   # pentacosan-12-ol
        "OCCCCCCCCCCCCO",             # 1,12-dodecanediol (should not be classified since it is diol, not monoalcohol)
    ]
    for smi in test_smiles:
        result, reason = is_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nClassified as fatty alcohol? {result}\nReason: {reason}\n{'-'*60}")