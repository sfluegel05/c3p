"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
Below is one improved solution. My reasoning was as follows. The previous code used a depth‐first search from the “alpha carbon” (the first carbon attached to a carboxyl group) without extra checks for the molecule’s “fatty acyl” character. That approach sometimes “found” a long continuous carbon path in very complex molecules (for example, in macrolides or polycyclic compounds with side‐chains) and called them fatty acids even though they were not simple linear acids. To better “capture” an ultra‐long–chain fatty acid (that is, a molecule whose carbons derive almost entirely from a single unbranched or nearly unbranched alkyl chain attached to –COOH) I made two changes. First, during the DFS I restrict the search to carbons that are not in rings (so that cyclic portions or fused polycyclic parts are left out) and—as before only follow carbon–carbon bonds. Second, after finding the longest chain starting from the alpha carbon I also compare that count to the total number of carbons in the molecule. In a true fatty acid the vast majority of carbons should be part of that chain (aside from minor substituents such as a hydroxyl substituent or a methyl branch). If extra carbon “limbs” or fragments are present, the molecule is likely a complex structure (and a false positive). In my code I require that the longest chain (plus the carboxyl carbon) comes within three carbons of the total carbon count. (This threshold is adjustable.) Finally, the threshold for “ultra‐long–chain” is that the total chain length must be greater than 27 (i.e. at least 28 carbons).

Below is the complete code in a single Python block:

------------------------------------------------

#!/usr/bin/env python
"""
Classifies: ultra-long-chain fatty acid
Definition: Any very long-chain fatty acid is defined as one where the fatty acid portion,
i.e. the carbon chain attached to a carboxylic acid (-COOH) group, consists largely of a single,
continuous (acyclic) chain whose total number of carbons (including the carboxyl carbon)
is greater than 27. A true fatty acid typically has only minimal additional substituents.
This classifier therefore (i) confirms the presence of a carboxylic acid group,
(ii) identifies the longest continuous non-ring carbon chain starting from the alpha carbon attached to the carboxyl,
and (iii) checks that nearly all carbons in the molecule are part of that chain.
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid, defined as having a total (acid) chain length > C27.
    
    The strategy is as follows:
      1. Parse the SMILES string.
      2. Identify a carboxylic acid group using a SMARTS pattern.
      3. From the carboxyl carbon (the carbonyl carbon) identify an adjacent carbon (the alpha carbon).
      4. Traverse only through non-ring carbon atoms starting with the alpha carbon to obtain the longest continuous chain.
      5. Include the carboxyl carbon in the count.
      6. To better ensure that the molecule is a simple fatty acid (and not a complex molecule with extra carbon fragments),
         compare the total number of carbons in the molecule to the length of the detected chain.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule is classified as an ultra-long-chain fatty acid, False otherwise.
      str: The reason for classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS to match a carboxylic acid group: carbonyl C(=O) with an -OH.
    carboxyl_smarts = "[CX3](=O)[OX1H]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:
        return False, "Molecule does not contain a carboxylic acid group"
    
    # Use the first match.
    # In our SMARTS pattern, we assume atom 0 is the carbonyl carbon.
    carboxyl_idx = matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Look for an alpha carbon attached to the carboxyl carbon.
    alpha_idx = None
    for neighbor in carboxyl_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # carbon
            alpha_idx = neighbor.GetIdx()
            break
    if alpha_idx is None:
        return False, "Carboxyl carbon is not attached to any carbon (no alkyl chain found)"
    
    # We now want to determine the longest continuous chain starting from the alpha carbon.
    # For a fatty acid, the chain should be largely acyclic. Thus, restrict to carbon atoms
    # that are not in rings.
    #
    # We'll use a DFS (depth-first-search) that only follows neighbors that are carbon and not in a ring.
    def dfs(atom_idx, visited):
        max_length = 1  # count the current carbon
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
    # so we do not backtrack into the carboxyl group.
    chain_length_from_alpha = dfs(alpha_idx, {carboxyl_idx, alpha_idx})
    # Total chain length: include the carboxyl carbon.
    total_chain_length = 1 + chain_length_from_alpha

    # For an ultra-long-chain fatty acid, the total chain length must be > 27.
    if total_chain_length <= 27:
        return False, f"Chain length is {total_chain_length} carbons, which is not greater than C27"
    
    # Now, ensure that the fatty acid chain makes up nearly the entire carbon content of the molecule.
    # Count the total number of carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Allow for minor substituents (such as a methyl branch or hydroxyl groups) by tolerating a difference of up to 3 carbons.
    if total_carbons - total_chain_length > 3:
        return False, (f"Chain length is {total_chain_length} carbons but the molecule has {total_carbons} carbons. "
                       "Excess carbon fragments suggest it is not a simple fatty acid.")
    
    return True, f"Chain length is {total_chain_length} carbons, qualifies as ultra-long-chain fatty acid"

# Example testing (uncomment the lines below to run a test)
# if __name__ == "__main__":
#     test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"  # example: dotriacontanoic acid
#     result, reason = is_ultra_long_chain_fatty_acid(test_smiles)
#     print(result, reason)
