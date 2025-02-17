"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: long‐chain fatty acid (chain length from C13 to C22)

A long‐chain fatty acid is defined here as a molecule that has a carboxylic acid group 
whose “acid carbon” is connected to a main (i.e. longest) carbon‐chain with a total count 
of between 13 and 22 carbons (inclusive). In our heuristic the acid carbon must be attached 
to one carbon. We then search the carbon–only graph for the longest simple path starting from 
the acid carbon. Finally, along that path we allow extra carbon substituents only if they themselves 
are terminal methyl groups (degree = 1).

Note: This is a heuristic method that will not resolve every “borderline” case.
"""

from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long‐chain fatty acid based on its SMILES string.
    Procedure:
      1. Look for one or more carboxylic acid groups (either –COOH or –COO–).
      2. For each matching acid group, identify the acid carbon (defined here as the C in C(=O)[OX2H,OX1-]).
      3. Require that the acid carbon has exactly one carbon neighbor (the start of the chain).
      4. Use a depth‐first search (DFS) on the sub‐graph of carbon atoms (ignoring all else) to find 
         the longest simple (acyclic) path starting from the acid carbon.
      5. For every atom in that candidate chain, verify that any additional carbon substituent that is 
         not part of the chosen chain is “small” (i.e. a terminal methyl group with degree 1). 
      6. If one acid group gives a chain with total carbon count (acid carbon plus chain) between 13 
         and 22 (inclusive) with no disallowed branching, classify as long‐chain fatty acid.
    
    Args:
      smiles (str): SMILES string of the molecule

    Returns:
      bool: True if the molecule satisfies our criteria, False otherwise
      str: A reason explaining the classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS for a carboxylic acid group: acid carbon is the one double-bonded to O and single-bonded to O (which may be protonated or anionic).
    acid_smarts = "[CX3](=O)[OX2H,OX1-]"
    acid_group = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_group)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # DFS function to find the longest simple path (list of atom indices) through carbon atoms
    def dfs_longest_path(current_idx, visited):
        current_atom = mol.GetAtomWithIdx(current_idx)
        best_path = [current_idx]
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:  # only consider carbons
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            new_visited = visited | {nbr_idx}
            sub_path = dfs_longest_path(nbr_idx, new_visited)
            candidate = [current_idx] + sub_path
            if len(candidate) > len(best_path):
                best_path = candidate
        return best_path

    messages = []
    for match in acid_matches:
        acid_idx = match[0]  # in our SMARTS, first atom is the acid carbon
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        # Acid carbon should have exactly one carbon neighbor (the beginning of the chain)
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            messages.append("Acid carbon does not have exactly one carbon neighbor (found {})".format(len(carbon_neighbors)))
            continue
        chain_start = carbon_neighbors[0]
        # Get the longest carbon chain starting from the acid carbon.
        # Include the acid carbon in the chain.
        longest_chain = dfs_longest_path(acid_idx, {acid_idx})
        chain_length = len(longest_chain)
        
        # Now check for disallowed branching.
        # We allow a substituent as long as it is a carbon atom that is terminal (degree 1).
        branching_disqualified = False
        for atom_idx in longest_chain:
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() != 6:
                    # ignore non-carbon neighbours
                    continue
                nbr_idx = nbr.GetIdx()
                # if the neighbor is not in the longest chain then it is a substituent
                if nbr_idx not in longest_chain:
                    # Allow if this branch is just one carbon (terminal methyl) 
                    if nbr.GetDegree() == 1:
                        continue
                    else:
                        branching_disqualified = True
                        break
            if branching_disqualified:
                break

        if branching_disqualified:
            messages.append("Branching detected along the candidate chain (non-terminal carbon substituent found)")
            continue

        # Check that the chain length falls in our accepted range: note that we count the acid carbon too.
        if chain_length < 13:
            messages.append("Longest carbon chain is {} carbons, which is too short".format(chain_length))
            continue
        if chain_length > 22:
            messages.append("Longest carbon chain is {} carbons, which is too long".format(chain_length))
            continue

        return True, "Contains a carboxylic acid group with a connected carbon chain of {} carbons".format(chain_length)
    
    if messages:
        return False, messages[0]
    else:
        return False, "No suitable carboxylic acid group with a connected carbon chain found"

# Example usage (uncomment to test):
# test_smiles = "OC(=O)CCCCCCCCCCCCCCCF"  # 16-Fluorohexadecanoic acid: should be True
# result, reason = is_long_chain_fatty_acid(test_smiles)
# print(result, reason)