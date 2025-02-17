"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: ultra-long-chain fatty acid
Definition: A very long-chain fatty acid is defined as one where the fatty acid portion --
that is, the carbon chain attached to a carboxylic acid (-COOH) group --
consists largely of a single, continuous (acyclic) chain whose total number of carbons
(including the carboxyl carbon) is greater than 27. A true fatty acid typically has only
minimal additional substituents. This classifier therefore:
  (i) confirms the presence of a carboxylic acid group,
  (ii) identifies the longest continuous non-ring carbon chain starting from the alpha carbon attached
      to the carboxyl, and
  (iii) checks that nearly all carbons in the molecule are part of that chain.
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid, defined as having a total (acid)
    chain length greater than C27.
    
    The strategy is as follows:
      1. Parse the SMILES string.
      2. Identify a carboxylic acid group using a SMARTS pattern.
      3. From the carboxyl carbon (the carbonyl carbon) identify an adjacent carbon (the alpha carbon).
      4. Traverse only through non-ring carbon atoms starting with the alpha carbon to obtain the longest continuous chain.
      5. Include the carboxyl carbon in the count.
      6. To ensure that the fatty acid chain makes up nearly the entire carbon content of the molecule,
         compare the total number of carbons in the molecule to the length of the detected chain.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as an ultra-long-chain fatty acid, False otherwise.
      str: The reason for classification.
    """
    
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS to match a carboxylic acid group: carbonyl C(=O) with a hydroxyl [OX1H]
    carboxyl_smarts = "[CX3](=O)[OX1H]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:
        return False, "Molecule does not contain a carboxylic acid group"
    
    # Use the first match; we assume atom 0 in the match is the carbonyl carbon.
    carboxyl_idx = matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Look for an alpha carbon attached to the carboxyl carbon (should be a carbon atom).
    alpha_idx = None
    for neighbor in carboxyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # carbon
            alpha_idx = neighbor.GetIdx()
            break
    if alpha_idx is None:
        return False, "Carboxyl carbon is not attached to any carbon (no alkyl chain found)"
    
    # Depth-First Search (DFS) to find the longest continuous chain starting from the alpha carbon.
    # We restrict the search to carbon atoms that are not in rings (acyclic).
    def dfs(atom_idx, visited):
        max_length = 1  # count the current carbon atom
        current_atom = mol.GetAtomWithIdx(atom_idx)
        for nbr in current_atom.GetNeighbors():
            # Only follow if neighbor is carbon, not visited, and not in a ring.
            if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() not in visited) and (not nbr.IsInRing()):
                new_visited = visited | {nbr.GetIdx()}
                branch_length = 1 + dfs(nbr.GetIdx(), new_visited)
                if branch_length > max_length:
                    max_length = branch_length
        return max_length

    # Begin DFS from the alpha carbon. Mark both the carboxyl and alpha carbons as visited
    # so as not to backtrack into the carboxyl group.
    chain_length_from_alpha = dfs(alpha_idx, {carboxyl_idx, alpha_idx})
    # Total chain length: include the carboxyl carbon.
    total_chain_length = 1 + chain_length_from_alpha

    # For an ultra-long-chain fatty acid, the total chain length must be greater than 27.
    if total_chain_length <= 27:
        return False, f"Chain length is {total_chain_length} carbons, which is not greater than C27"
    
    # Count the total number of carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Allow for minor substituents by tolerating a difference of up to 3 carbons.
    if total_carbons - total_chain_length > 3:
        return False, (f"Chain length is {total_chain_length} carbons but the molecule has {total_carbons} carbons. "
                       "Excess carbon fragments suggest it is not a simple fatty acid.")
    
    return True, f"Chain length is {total_chain_length} carbons, qualifies as ultra-long-chain fatty acid"

# Example testing (uncomment to run a test)
# if __name__ == "__main__":
#     test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"  # example: dotriacontanoic acid
#     result, reason = is_ultra_long_chain_fatty_acid(test_smiles)
#     print(result, reason)