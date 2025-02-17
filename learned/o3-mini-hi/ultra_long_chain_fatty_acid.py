"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: ultra-long-chain fatty acid
Definition: Any very long-chain fatty acid which has a chain length greater than C27.
A fatty acid is identified as a molecule containing a carboxyl group (–COOH)
with an attached linear alkyl chain.
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid, defined as having a total chain length > C27.
    The chain is defined as the continuous chain of carbons connected to the carboxylic acid group.
    
    For our purpose, we require that the molecule contains a carboxylic acid substructure,
    and we then search for the longest continuous carbon chain (starting from the alpha carbon attached
    to the carboxyl carbon) using a depth-first search. The chain length (number of carbons in the chain,
    including the carboxyl carbon) is then compared against the threshold of 27.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a SMARTS pattern for a carboxylic acid group.
    # This pattern matches a carbon (index 0) double bonded to an oxygen and single bonded to an -OH.
    ca_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(ca_pattern)
    if not matches:
        return False, "Molecule does not contain a carboxylic acid group"
    
    # Use the first occurrence of the carboxyl group.
    # In our SMARTS "C(=O)O", atom index 0 is the carbonyl carbon.
    carboxyl_idx = matches[0][0]
    
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    # Look for a neighboring carbon (alpha carbon) attached to the carboxyl carbon.
    alpha_idx = None
    for neighbor in carboxyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # carbon
            alpha_idx = neighbor.GetIdx()
            break
    if alpha_idx is None:
        return False, "Carboxyl carbon is not attached to any carbon (no alkyl chain found)"
    
    # We now perform a DFS on the molecule's carbon-only substructure.
    # Our goal is to find the length (number of carbon atoms) of the longest continuous chain that originates
    # from the alpha carbon. We then add the carboxyl carbon (from which the chain starts) to get total acid chain length.
    def dfs(atom_idx, visited):
        # Returns the maximum chain length (number of carbons in this branch including the current atom)
        max_length = 1
        current_atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in current_atom.GetNeighbors():
            # Only follow if neighbor is carbon and not yet visited.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in visited:
                new_visited = visited | {nbr.GetIdx()}
                branch_length = 1 + dfs(nbr.GetIdx(), new_visited)
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    # Start DFS from the alpha carbon. We start visited with both the carboxyl carbon and the alpha carbon
    # so the search does not go back to the carboxyl group.
    chain_length_from_alpha = dfs(alpha_idx, {carboxyl_idx, alpha_idx})
    # Total chain length includes the carboxyl carbon.
    total_chain_length = 1 + chain_length_from_alpha

    if total_chain_length <= 27:
        return False, f"Chain length is {total_chain_length} carbons, which is not greater than C27"
    
    return True, f"Chain length is {total_chain_length} carbons, qualifies as ultra-long-chain fatty acid"

# Example testing (uncomment to run)
# if __name__ == "__main__":
#     # Test with one of the provided examples (dotriacontanoic acid)
#     test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"
#     result, reason = is_ultra_long_chain_fatty_acid(test_smiles)
#     print(result, reason)