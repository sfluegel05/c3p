"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: ultra-long-chain fatty acid
Definition: A very long-chain fatty acid is defined as one where the fatty acid portion --
that is, the carbon chain attached to the carboxylic acid (-COOH) group --
consists largely of a single, continuous (acyclic) chain whose total number of carbons
(including the carboxyl carbon) is greater than 27. A true fatty acid typically has only
minimal additional substituents.
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid, defined as having a total (acid)
    chain length greater than C27.

    The strategy is:
      1. Parse the SMILES string.
      2. Identify a carboxylic acid group using an updated SMARTS pattern.
      3. From the carboxyl carbon (the carbonyl carbon) identify an adjacent carbon atom (alpha carbon).
      4. Perform a Depth-First Search (DFS) on non-ring carbons from the alpha carbon to find the longest continuous chain.
      5. Include the carboxyl carbon in the overall fatty acid chain length.
      6. Check if nearly all carbons in the molecule belong to this chain (allowing a small difference).

    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule qualifies as an ultra-long-chain fatty acid, False otherwise.
      str: Explanation for the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Updated SMARTS for a carboxylic acid group:
    # The pattern [CX3](=O)[OX2H1] matches a trivalent carbon (as in a carbonyl) double-bonded to O 
    # and single-bonded to an -OH group, where the hydroxyl oxygen is represented as [OX2H1].
    carboxyl_smarts = "[CX3](=O)[OX2H1]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:
        return False, "Molecule does not contain a carboxylic acid group"
    
    # Assume the first match; the first atom in the match is taken as the carboxyl carbon.
    carboxyl_idx = matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Look for an alpha carbon attached to the carboxyl carbon. It must be a carbon atom.
    alpha_idx = None
    for neighbor in carboxyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Only process carbon atoms.
            alpha_idx = neighbor.GetIdx()
            break
    if alpha_idx is None:
        return False, "Carboxyl carbon is not attached to any carbon (no alkyl chain found)"
    
    # Recursive Depth-First Search to find the longest continuous chain from a starting carbon atom.
    # We restrict the search to non-ring (acyclic) carbon atoms.
    def dfs(atom_idx, visited):
        max_length = 1  # Count the current carbon
        current_atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in current_atom.GetNeighbors():
            if (nbr.GetAtomicNum() == 6 and (nbr.GetIdx() not in visited)
                and (not nbr.IsInRing())):
                new_visited = visited | {nbr.GetIdx()}
                branch_length = 1 + dfs(nbr.GetIdx(), new_visited)
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    # Begin DFS from the alpha carbon.
    # Mark both the carboxyl and alpha carbons as visited to avoid backtracking.
    chain_length_from_alpha = dfs(alpha_idx, {carboxyl_idx, alpha_idx})
    # Total chain length includes the carboxyl carbon.
    total_chain_length = 1 + chain_length_from_alpha

    # Check if the chain is long enough: strictly greater than 27 carbons.
    if total_chain_length <= 27:
        return False, f"Chain length is {total_chain_length} carbons, which is not greater than C27"
    
    # Count the total number of carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Allow a difference of up to 3 carbons between the chain length and the total carbon count.
    if total_carbons - total_chain_length > 3:
        return False, (f"Chain length is {total_chain_length} carbons but the molecule has {total_carbons} carbons. "
                       "Excess carbon fragments suggest it is not a simple fatty acid.")
    
    return True, f"Chain length is {total_chain_length} carbons, qualifies as ultra-long-chain fatty acid"

# Example testing (uncomment to run a test)
# if __name__ == "__main__":
#     test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"  # Example: dotriacontanoic acid
#     result, reason = is_ultra_long_chain_fatty_acid(test_smiles)
#     print(result, reason)