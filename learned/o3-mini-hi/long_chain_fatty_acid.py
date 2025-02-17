"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: long-chain fatty acid (chain length from C13 to C22)

A long-chain fatty acid is defined here as a molecule that has a carboxylic acid group 
whose “acid carbon” is connected to a main (i.e. longest) carbon‐chain with between 13 
and 22 carbons. In our heuristic the acid carbon must have exactly one carbon neighbor, 
and the candidate chain is grown on the carbon-only graph (excluding ring atoms and 
disallowing any branching beyond a terminal methyl substituent). In addition, we require 
that the candidate chain be the dominant carbon substructure (covering at least 80% of 
the carbons in the connected component).
  
Note: This is a heuristic and may still fail for borderline or highly complex cases.
"""

from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long‐chain fatty acid based on its SMILES string.
    
    Strategy:
      1. Identify carboxylic acid groups using SMARTS "[CX3](=O)[OX2H,OX1-]". The acid carbon (first atom)
         is expected to be connected to exactly one carbon neighbor (the beginning of the chain).
      2. Starting from that neighbor (with the acid carbon included) use DFS on the carbon-only graph 
         (only atoms with atomic number 6) to enumerate all simple paths. 
      3. Reject any candidate chain if any carbon in the chain has another carbon neighbor outside the 
         chain (i.e. branching) or if any atom in the candidate lies in a ring.
      4. Among the valid candidate chains, pick the longest one. Then, if the chain has between 13 and 22 
         carbons (inclusive), further require that it “dominates” the connected carbon subgraph (i.e. it 
         accounts for at least 80% of the carbons in that fragment). If so, classify as a long‐chain fatty acid.
    
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      bool: True if molecule meets criteria, False otherwise
      str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS to find a carboxylic acid group; assuming the acid carbon is the first atom.
    acid_smarts = "[CX3](=O)[OX2H,OX1-]"
    acid_group = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_group)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Precalculate carbon-only neighbors for every atom (only count neighbors that are carbon)
    carbon_neighbors_dict = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            carbon_neighbors_dict[atom.GetIdx()] = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    
    # Helper: get connected component of carbons (the set of carbon indices) starting from a given index.
    def get_carbon_component(start_idx):
        to_visit = [start_idx]
        comp = set()
        while to_visit:
            idx = to_visit.pop()
            if idx in comp:
                continue
            comp.add(idx)
            for nbr in carbon_neighbors_dict.get(idx, []):
                if nbr not in comp:
                    to_visit.append(nbr)
        return comp
    
    # DFS generator to yield all simple paths (list of atom indices) in the carbon-only graph
    def dfs_paths(current_idx, path, visited):
        # Always yield the current path so far.
        yield path
        for nbr in carbon_neighbors_dict.get(current_idx, []):
            if nbr in visited:  # avoid cycles
                continue
            new_visited = visited | {nbr}
            yield from dfs_paths(nbr, path + [nbr], new_visited)
    
    messages = []
    # Iterate over each carboxylic acid match.
    for match in acid_matches:
        acid_idx = match[0]  # acid carbon as per our SMARTS
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        # The acid carbon must be connected to exactly one carbon (the start of the fatty acid chain)
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            messages.append("Acid carbon does not have exactly one carbon neighbor (found {})".format(len(carbon_neighbors)))
            continue
        chain_start = carbon_neighbors[0]
        # Start candidate chain with the acid carbon and its unique carbon neighbor.
        initial_path = [acid_idx, chain_start]
        
        # Get all candidate chains starting at chain_start (paths always include acid_idx at beginning).
        candidate_paths = list(dfs_paths(chain_start, initial_path, set(initial_path)))
        if not candidate_paths:
            messages.append("No candidate chain found starting from the carboxylic acid group")
            continue
        
        # Filter candidate chains: ensure no branch (if any extra carbon neighbor not in the chain is found)
        # and ensure none of the atoms in the chain are in a ring.
        valid_candidates = []
        for path in candidate_paths:
            branch_found = False
            for idx in path:
                atom = mol.GetAtomWithIdx(idx)
                # Reject chain if any atom is in a ring.
                if atom.IsInRing():
                    branch_found = True
                    break
                # Check each carbon neighbor – if it is not in the current chain, then it is a branch.
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() != 6:
                        continue
                    nbr_idx = nbr.GetIdx()
                    if nbr_idx not in path:
                        # Allow if this neighbor is terminal (only one carbon neighbor overall)
                        if len(carbon_neighbors_dict.get(nbr_idx, [])) != 1:
                            branch_found = True
                            break
                if branch_found:
                    break
            if not branch_found:
                valid_candidates.append(path)
        
        if not valid_candidates:
            messages.append("Branching (or ring-containing atoms) detected along the candidate chain")
            continue
        
        # Pick best candidate chain (the longest)
        best_path = max(valid_candidates, key=lambda p: len(p))
        chain_length = len(best_path)
        if chain_length < 13:
            messages.append("Longest carbon chain is {} carbons, which is too short".format(chain_length))
            continue
        if chain_length > 22:
            messages.append("Longest carbon chain is {} carbons, which is too long".format(chain_length))
            continue
        
        # Ensure the candidate chain is the dominant carbon substructure.
        # Get the connected component (in the carbon-only graph) containing the acid carbon.
        comp = get_carbon_component(acid_idx)
        ratio = chain_length / len(comp)  if comp else 0
        if ratio < 0.8:
            messages.append("Candidate chain ({} carbons) makes up only {:.0%} of the connected carbon component".format(chain_length, ratio))
            continue
        
        return True, "Correct: Contains a carboxylic acid group with a connected carbon chain of {} carbons".format(chain_length)
    
    # If no acid group yielded a valid candidate, return the first message.
    if messages:
        return False, messages[0]
    else:
        return False, "No suitable carboxylic acid group with an attached carbon chain found"

# Example usage (uncomment to test):
# test_smiles = "OC(=O)CCCCCCCCCCCCCCCF"  # 16-Fluorohexadecanoic acid (should be True)
# result, reason = is_long_chain_fatty_acid(test_smiles)
# print(result, reason)