"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: ultra-long-chain fatty acid
Definition: A very long-chain fatty acid is defined as one where the fatty acid portion – that is, the carbon chain attached to the carboxylic acid (-COOH) group – 
consists largely of a single, continuous (acyclic and mostly unbranched) chain whose total number of carbons (including the carboxyl carbon) is greater than 27.
We further require that the molecule contains exactly one carboxylic acid group, that the carboxyl carbon is terminal (attached to only a single alkyl carbon),
and that aside from an optionally allowed hydroxyl at the alpha carbon the main chain is not substituted with extra heteroatoms.
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid.
    
    Strategy:
      1. Parse the SMILES and check that exactly one carboxylic acid group (SMARTS "[CX3](=O)[OX2H1]") is present.
      2. Ensure that the acid carbon (the first atom in the match) is terminal – meaning it has only one carbon neighbor.
      3. From that neighbor (the alpha carbon), perform a DFS across carbon atoms (only considering acyclic atoms) 
         to find the longest simple (non-repeating) carbon path. (Branching is allowed in the search but not in the final “main chain”.)
      4. Compute the total chain length = (1 for carboxyl carbon + length of longest carbon path).
         The chain must have >27 carbons.
      5. Check that almost all carbon atoms in the molecule are on that main chain (no extensive extra fragments).
      6. Finally, check that the main chain is “clean”: aside from a possible hydroxyl substitution on the alpha carbon,
         no other chain carbons are substituted by oxygen (for example, the terminal carbon must be a plain CH3).
         
    Args:
      smiles (str): SMILES string for the molecule.
      
    Returns:
      bool: True if the molecule qualifies as an ultra-long-chain fatty acid, otherwise False.
      str: Explanation for the classification decision.
      
    If the task cannot be done unambiguously the function may return (None, None).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1. Find carboxylic acid groups using SMARTS.
    carboxyl_smarts = "[CX3](=O)[OX2H1]"
    carboxyl_query = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol.GetSubstructMatches(carboxyl_query)
    if len(carboxyl_matches) != 1:
        return False, f"Expected exactly one carboxylic acid group; found {len(carboxyl_matches)}"

    # Assume first (and only) match:
    # In the SMARTS, the first atom is the carboxyl (carbonyl) carbon.
    carboxyl_idx = carboxyl_matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)

    # Step 2. Check that the carboxyl carbon is terminal (should have exactly one carbon neighbor)
    carbon_neighbors = [nbr for nbr in carboxyl_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "Carboxyl carbon is not terminal (expected one carbon neighbor)"
    alpha_atom = carbon_neighbors[0]
    alpha_idx = alpha_atom.GetIdx()

    # Step 3. Traverse the carbon chain (only atomic number 6) starting from the alpha carbon.
    # We do a DFS that finds the longest simple path. We allow a branch but we return only the longest chain.
    def dfs(current_idx, visited):
        current_atom = mol.GetAtomWithIdx(current_idx)
        max_length = 1
        max_path = [current_idx]
        # Consider neighbors that are carbons – we allow unsaturation.
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            n_idx = nbr.GetIdx()
            if n_idx in visited:
                continue
            # Only consider if the atom is acyclic (we skip those in rings) to be sure we get a fatty acid chain
            if nbr.IsInRing():
                continue
            new_visited = visited.copy()
            new_visited.add(n_idx)
            rec_length, rec_path = dfs(n_idx, new_visited)
            if 1 + rec_length > max_length:
                max_length = 1 + rec_length
                max_path = [current_idx] + rec_path
        return max_length, max_path

    # Start DFS from the alpha carbon; include the carboxyl carbon in visited so we do not go back.
    visited = set([carboxyl_idx, alpha_idx])
    chain_length_from_alpha, best_path = dfs(alpha_idx, visited)
    # Our total chain includes the carboxyl carbon at the head.
    total_chain_length = 1 + chain_length_from_alpha  # carboxyl carbon + chain
    
    if total_chain_length <= 27:
        return False, f"Chain length is {total_chain_length} carbons, which is not greater than C27"
    
    # The main chain we propose is: [carboxyl_atom] + best_path (which begins with alpha)
    main_chain = [carboxyl_idx] + best_path

    # Step 4. Check that most (or all) carbon atoms belong to this chain.
    # Allow at most one extra carbon (to allow e.g. a terminal methyl branch in a genuine fatty acid)
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons - len(main_chain) > 1:
        return False, (f"Chain length is {len(main_chain)} carbons but molecule has {total_carbons} carbons. "
                       "Extra carbon fragments suggest it is not a simple fatty acid.")
    
    # Step 5. Examine the main chain for extra (non-allowed) substitutions.
    # Allowed: the carboxyl carbon (head) and the alpha carbon may have a hydroxyl substituent.
    # All other carbons in the main chain should not have any additional oxygen neighbors.
    # Also the terminal carbon must be a plain CH3 (i.e. no oxygen attached).
    for pos, atom_idx in enumerate(main_chain):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Gather oxygen neighbors that are not in the main chain.
        oxy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in main_chain]
        # Allow oxygen at the alpha carbon (position 1 in main_chain after carboxyl) only.
        if pos == 0:
            # carboxyl carbon: already engaged in acid group
            continue
        elif pos == 1:
            # alpha carbon: allowed to have one OH substituent.
            if len(oxy_neighbors) > 1:
                return False, f"Alpha carbon has too many oxygen substituents"
        else:
            # Other positions: no extra oxygen allowed.
            if oxy_neighbors:
                return False, f"Extra oxygen substituent found at chain position {pos+1}"
    # Check that the terminal carbon (last in main_chain) is a CH3 (i.e. no oxygen neighbors)
    terminal_atom = mol.GetAtomWithIdx(main_chain[-1])
    term_oxy = [nbr for nbr in terminal_atom.GetNeighbors() if nbr.GetAtomicNum()==8]
    if term_oxy:
        return False, "Terminal carbon is substituted with oxygen (should be a terminal CH3)"
    
    # If all checks pass, qualify as an ultra-long-chain fatty acid.
    return True, f"Chain length is {total_chain_length} carbons, qualifies as ultra-long-chain fatty acid"

# Example usage (if run as script, uncomment below):
# if __name__=="__main__":
#     examples = [
#         "CCCCCC\\C=C/CCCCCCCCCC(O)C(O)=O",  # Example: (12Z,15Z,18Z,21Z,24Z)-2-hydroxytriacontapentaenoic acid-like
#         "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O",  # A simple unbranched fatty acid (check chain length)
#     ]
#     for smi in examples:
#         res, reason = is_ultra_long_chain_fatty_acid(smi)
#         print(smi, res, reason)