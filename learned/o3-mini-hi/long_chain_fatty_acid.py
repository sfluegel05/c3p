"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: long‐chain fatty acid (chain length from C13 to C22)
A long‐chain fatty acid is defined here as a molecule that has a carboxylic acid group 
whose “acid carbon” is connected to a main (i.e. longest) carbon‐chain with a total count 
of between 13 and 22 carbons (inclusive). In our heuristic we also require that no chain carbon
bears an extra carbon substituent (i.e. the main chain must be unbranched).
Note: This is a heuristic method that will not resolve every “borderline” case.
"""

from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long‐chain fatty acid based on its SMILES string.
    The procedure is:
     1. Find one or more carboxylic acid groups. We accept both protonated (–COOH) and 
        deprotonated (–COO–) forms.
     2. For each carboxyl match, identify the “acid carbon” (the one in C(=O)–O).
     3. Check that the acid carbon is attached to exactly one carbon atom.
     4. Do a depth‑first search (DFS) over the carbon sub‐graph to find the longest simple 
        (acyclic) path (this is taken as the “alkyl chain” of the acid).
     5. Finally, for every carbon in that chain, check that it does not have any extra carbon
        substituents (other than the one coming from the chain itself).
    If one acid group gives a main chain with total carbon count (acid carbon plus chain) between
    13 and 22 and no side–chain branching, then it returns True with an explanation.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule conforms to our long‐chain fatty acid criteria, otherwise False.
        str: A reason explaining the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS for a carboxylic acid group.
    # This pattern matches a carbon with a double-bonded oxygen and single-bonded oxygen (which may be protonated or anionic).
    acid_smarts = "[CX3](=O)[OX2H,OX1-]"
    acid_group = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_group)
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # Helper function: use DFS to find all simple paths (over carbon atoms) starting from a given atom.
    def dfs_longest_path(current_idx, visited):
        current_atom = mol.GetAtomWithIdx(current_idx)
        best_path = [current_idx]
        # Consider only neighbours that are carbons.
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # Extend the visited set and explore further.
            new_visited = visited | {nbr_idx}
            sub_path = dfs_longest_path(nbr_idx, new_visited)
            candidate = [current_idx] + sub_path
            if len(candidate) > len(best_path):
                best_path = candidate
        return best_path

    # For each carboxyl acid match, try to get a chain that meets our criteria.
    messages = []
    for match in acid_matches:
        # In our SMARTS the first atom is the acid carbon.
        acid_idx = match[0]
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        # Find carbon neighbours of the acid carbon.
        carbon_nei = [nbr.GetIdx() for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_nei) != 1:
            messages.append("Acid carbon does not have exactly one carbon neighbor (found {})".format(len(carbon_nei)))
            continue
        chain_start = carbon_nei[0]
        # Compute the longest carbon chain that can be reached from the acid carbon.
        # We include the acid carbon in the count.
        longest_chain = [acid_idx] + dfs_longest_path(chain_start, {acid_idx, chain_start})
        chain_length = len(longest_chain)
        # Now check that none of the chain carbons have an extra carbon substituent.
        # For each atom in the chain except (possibly) the terminal one, check that among its carbon neighbours,
        # the ones outside the chain are null.
        branching_flag = False
        for idx in longest_chain:
            atom = mol.GetAtomWithIdx(idx)
            # Get all carbon neighbours
            c_neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            # Allowed neighbours are those already in the chain.
            for nbr_idx in c_neighbors:
                if nbr_idx not in longest_chain:
                    branching_flag = True
                    break
            if branching_flag:
                break

        if branching_flag:
            messages.append("Branching detected along the chain starting at acid carbon")
            continue

        if chain_length < 13:
            messages.append("Longest carbon chain is {} carbons, which is too short".format(chain_length))
            continue
        if chain_length > 22:
            messages.append("Longest carbon chain is {} carbons, which is too long".format(chain_length))
            continue

        # Passed all checks: we have an acid group and an unbranched carbon chain of acceptable length.
        return (True, "Contains a carboxylic acid group with a connected unbranched carbon chain of {} carbons".format(chain_length))

    # If we got here then no acid group produced an acceptable chain.
    if messages:
        return (False, messages[0])
    return (False, "No suitable carboxylic acid group with a connected carbon chain found")

# Example (uncomment to test):
# test_smiles = "OC(=O)CCCCCCCCCCCCCCC([2H])([2H])[2H]"  # Palmitic acid(d3): 16 carbons
# result, reason = is_long_chain_fatty_acid(test_smiles)
# print(result, reason)