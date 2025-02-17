"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: long‐chain fatty acid (chain length from C13 to C22)

A long‐chain fatty acid is defined here as a molecule that has a carboxylic acid group 
whose “acid carbon” is connected to a main (i.e. longest) carbon‐chain with a total count 
of between 13 and 22 carbons (inclusive). In our heuristic the acid carbon must have exactly 
one carbon neighbor (the start of the chain). Then we search the carbon–only graph for 
all possible simple paths starting at the acid carbon and choose those for which any extra 
carbon substituent (found off the candidate chain) is “small” (i.e. a terminal methyl group). 
If one candidate yields a chain of acceptable length, the molecule is classified as a long–chain fatty acid.

Note: This is a heuristic method and may miss borderline or branched cases.
"""

from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long‐chain fatty acid based on its SMILES string.
    
    Strategy:
      1. Find carboxylic acid groups (SMARTS: [CX3](=O)[OX2H,OX1-]). We consider the acid carbon (first atom).
      2. The acid carbon must have exactly one carbon neighbor (that becomes the start of the main chain).
      3. From that starting point (including the acid carbon) we enumerate all possible simple paths 
         in the carbon-only graph (all atoms with atomic number 6) using DFS.
      4. For each candidate chain (path), check that for every carbon in the chain, any extra
         carbon neighbor (not in the chain) is allowed only if it is a terminal (has a carbon degree of 1).
      5. If a candidate exists with total chain length between 13 and 22 (inclusive) carbons, return True.
      
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      bool: True if molecule passes the criteria, False otherwise.
      str: A reason explaining the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Pre-calculate carbon-only degree for all atoms (only counting neighbors that are carbons).
    carbon_degrees = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            # Count only neighbors that are carbons.
            carbon_degrees[atom.GetIdx()] = sum(1 for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6)
    
    # Define SMARTS for a carboxylic acid group.
    acid_smarts = "[CX3](=O)[OX2H,OX1-]"
    acid_group = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_group)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # This list will store messages from each acid group attempt.
    messages = []
    
    # Helper: Generator to yield all simple paths (as lists of atom indices) in the carbon-only graph.
    def dfs_paths(current_idx, path, visited):
        yield path
        current_atom = mol.GetAtomWithIdx(current_idx)
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            new_visited = visited | {nbr_idx}
            yield from dfs_paths(nbr_idx, path + [nbr_idx], new_visited)
    
    # For each acid group match, try the DFS.
    for match in acid_matches:
        acid_idx = match[0]  # as per our SMARTS, first atom is the acid carbon
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        # The acid carbon must be connected to exactly one carbon (the chain starter).
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            messages.append("Acid carbon does not have exactly one carbon neighbor (found {})".format(len(carbon_neighbors)))
            continue
        chain_start = carbon_neighbors[0]
        # Build candidate chain starting from the acid carbon.
        # We start with a minimal chain including the acid carbon and its unique neighbor.
        initial_path = [acid_idx, chain_start]
        
        # Collect all candidate paths starting from acid_idx (paths will always begin with acid_idx).
        candidate_paths = list(dfs_paths(chain_start, initial_path, set(initial_path)))
        if not candidate_paths:
            messages.append("No candidate chain found starting from the carboxylic acid group")
            continue
        
        # Filter candidate paths: check for disallowed branching.
        valid_candidate_paths = []
        for path in candidate_paths:
            branch_problem = False
            for atom_idx in path:
                atom = mol.GetAtomWithIdx(atom_idx)
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() != 6:
                        continue
                    nbr_idx = nbr.GetIdx()
                    # If neighbor not in the chosen main chain, it’s a substituent.
                    if nbr_idx not in path:
                        # Allow only if that neighbor is terminal (in the carbon-only graph; degree == 1).
                        if carbon_degrees.get(nbr_idx, 0) != 1:
                            branch_problem = True
                            break
                if branch_problem:
                    break
            if not branch_problem:
                valid_candidate_paths.append(path)
                
        if not valid_candidate_paths:
            messages.append("Branching detected along the candidate chain (non-terminal carbon substituent found)")
            continue
        
        # Now, among valid candidate chains, find the one with the maximum length.
        # (Note: the acid carbon is included.)
        best_path = max(valid_candidate_paths, key=lambda p: len(p))
        chain_length = len(best_path)
        # Check if the chain length is in the accepted range.
        if chain_length < 13:
            messages.append("Longest carbon chain is {} carbons, which is too short".format(chain_length))
            continue
        if chain_length > 22:
            messages.append("Longest carbon chain is {} carbons, which is too long".format(chain_length))
            continue
        
        # If here, we have a valid candidate.
        return True, "Contains a carboxylic acid group with a connected carbon chain of {} carbons".format(chain_length)
    
    # If none of the acid groups gave a valid candidate, return the first message.
    if messages:
        return False, messages[0]
    else:
        return False, "No suitable carboxylic acid group with a connected carbon chain found"

# Example usage (uncomment to test):
# test_smiles = "OC(=O)CCCCCCCCCCCCCCCF"  # 16-Fluorohexadecanoic acid (should be True)
# result, reason = is_long_chain_fatty_acid(test_smiles)
# print(result, reason)