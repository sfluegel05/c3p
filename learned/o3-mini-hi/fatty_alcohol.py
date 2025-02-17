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
  2. Searches for an alcohol group where a saturated (or stereochemically enriched) carbon is attached to –OH.
  3. For each candidate –OH the alcohol carbon is required to be acyclic (not found in a ring) and, aside
     from the -OH, only has carbon neighbors.
  4. A connectivity graph of non‐aromatic, non‐ring carbon atoms is constructed.
  5. A DFS is used to “grow” the longest simple chain (path) from the alcohol carbon.
  6. Along that chain, each carbon is checked to ensure it does not have any extra neighbors that are heteroatoms.
  7. If one candidate –OH gives a longest chain with at least MIN_CHAIN_LENGTH carbons, the molecule is classified as a fatty alcohol.
  
Note: Branching (extra carbon attachments) is allowed.
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
    
    # Define a SMARTS query for an alcohol group where a non‐ring carbon is bonded to a hydroxyl (-OH).
    # This captures [C;!R] (a non‐ring carbon) connected to [OX2H].
    alcohol_smarts = "[C;!R][OX2H]"
    alcohol_query = Chem.MolFromSmarts(alcohol_smarts)
    alcohol_matches = mol.GetSubstructMatches(alcohol_query)
    if not alcohol_matches:
        return False, "No proper alcohol (-OH) group found on a non‐ring carbon"
    
    # Build a set of carbon atoms that are non‐aromatic and not in a ring.
    carbon_nodes = set()
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6 and (not atom.GetIsAromatic()) and (not atom.IsInRing()):
            carbon_nodes.add(atom.GetIdx())
    
    # Build a connectivity graph (dictionary) among these carbon atoms.
    carbon_graph = {idx: [] for idx in carbon_nodes}
    for idx in carbon_nodes:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in carbon_nodes:
                carbon_graph[idx].append(nbr.GetIdx())
    
    # DFS helper: returns the longest simple (non-repeating) path starting from 'current' node.
    def dfs_longest(current, visited):
        best_path = [current]
        for nbr in carbon_graph.get(current, []):
            if nbr not in visited:
                path = dfs_longest(nbr, visited | {nbr})
                if len(path) + 1 > len(best_path):
                    best_path = [current] + path
        return best_path

    # Helper function to check the "purity" of the candidate chain.
    # For each carbon in the chain, any neighbor that is not part of the chain must be hydrogen or carbon.
    def chain_is_pure(chain):
        chain_set = set(chain)
        for cid in chain:
            atom = mol.GetAtomWithIdx(cid)
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in chain_set:
                    # Allow hydrogen (atomic num 1) or carbon (atomic num 6).
                    if nbr.GetAtomicNum() not in {1, 6}:
                        return False
        return True

    # Use a minimum chain length requirement.
    MIN_CHAIN_LENGTH = 7  # Counting the alcohol carbon itself.
    
    best_chain_length = 0
    candidate_found = False
    candidate_reason = ""
    
    # Iterate over candidate alcohol groups.
    for match in alcohol_matches:
        # Each match: first atom is the candidate alcohol carbon; second is the hydroxyl oxygen.
        alc_c_idx = match[0]
        alc_o_idx = match[1]
        alc_atom = mol.GetAtomWithIdx(alc_c_idx)
        
        # Skip candidate if the alcohol carbon is in a ring.
        if alc_atom.IsInRing():
            continue
        
        # Check: apart from the -OH, all neighbors of the alcohol carbon should be carbons.
        valid = True
        for nbr in alc_atom.GetNeighbors():
            if nbr.GetIdx() == alc_o_idx:
                continue
            if nbr.GetAtomicNum() != 6:
                valid = False
                break
        if not valid:
            continue
        
        # Ensure the alcohol carbon is in our aliphatic carbon graph.
        if alc_c_idx not in carbon_graph:
            continue
        
        # Compute the longest contiguous chain (simple path) starting from the alcohol carbon.
        longest_path = dfs_longest(alc_c_idx, {alc_c_idx})
        chain_length = len(longest_path)
        if chain_length > best_chain_length:
            best_chain_length = chain_length
        
        # Check if the candidate chain is "pure" (no attached heteroatoms other than H/C).
        if not chain_is_pure(longest_path):
            continue
        
        if chain_length >= MIN_CHAIN_LENGTH:
            candidate_found = True
            candidate_reason = (f"Molecule contains an alcohol group with a "
                                f"{chain_length}-carbon contiguous acyclic aliphatic chain (chain path: {longest_path}).")
            break  # Accept the first candidate meeting the criteria.

    if candidate_found:
        return True, candidate_reason
    else:
        return False, (f"Longest acyclic aliphatic chain from any candidate alcohol group is only "
                       f"{best_chain_length} carbons (need at least {MIN_CHAIN_LENGTH}).")
    
# Example usage:
if __name__ == "__main__":
    test_smiles = [
        "CCCCCCC(C)O",      # octan-2-ol (should be true)
        "CCCC[C@H](C)O",     # (2S)-2-heptanol
        "OC(CCCCCCCCCCCCCCC)CC(=O)C1=CC=CC=C1",  # 3-Hydroxy-1-phenyl-1-octadecanone
        "OC1=CC(O)=CC(=C1)CCCCCCCCCCCCCCC[C@@H](O)C",  # Hansfordiol A
        "O=C(O[C@@H](CCCCCCCCCCCCC[C@H](O)C)[C@H](C1=C)C(=O)O"  # Allo-murolic acid (with extra acid groups)
    ]
    for smi in test_smiles:
        result, reason = is_fatty_alcohol(smi)
        print(f"SMILES: {smi}\nClassified as fatty alcohol? {result}\nReason: {reason}\n{'-'*60}")