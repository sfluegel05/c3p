"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid (chain length > C22; > C27 = ultra-long-chain)
Heuristic: The molecule must have a carboxylic acid group. Then, starting from the acid carbon,
we identify its unique carbon neighbor and “walk” down the chain recursively.
At each step the chain must be unbranched (i.e. the current carbon, excluding the carbon we came from,
must have exactly 0 or 1 carbon neighbors). If we obtain a linear chain (including the acid carbon) longer than 22 carbons,
the molecule qualifies. We return a note if the chain is >27 carbons (ultra-long-chain).
"""

from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as having an unbranched fatty acid chain (starting
    from the carboxylic acid carbon) with more than 22 contiguous carbon atoms.
    
    The algorithm:
      1. Confirm that a carboxylic acid group (protonated or deprotonated) is present.
      2. Identify its acid carbon. Then from that acid carbon pick its unique carbon neighbor.
         If there is not exactly one such neighbor, we cannot clearly define a single fatty acid chain.
      3. Recursively “walk” down the chain. At each carbon (other than the one we came from),
         there must be exactly one other carbon to continue the linear (unbranched) chain.
         If 0 carbons are found, we have reached the terminus; if >1, then the chain is branched.
      4. Count the total chain length (including the acid carbon). If this length is >22,
         the molecule qualifies as a very long-chain fatty acid.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule qualifies as a very long-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find carboxylic acid group (allow protonated or deprotonated forms)
    acid_smarts1 = Chem.MolFromSmarts("C(=O)[OH]")
    acid_smarts2 = Chem.MolFromSmarts("C(=O)[O-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts1) + mol.GetSubstructMatches(acid_smarts2)
    
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Use the first match: assume the first atom is the acid carbon
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Among acid carbon's neighbors, select only carbon atoms.
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, f"Expected exactly one carbon neighbor off the acid carbon, found {len(carbon_neighbors)}"
    
    # Recursive function to walk along an unbranched chain.
    def walk_chain(atom, coming_from_idx):
        """
        Walks from the current atom along the chain.
        Returns the length of the chain from this atom (counting this atom as 1)
        if the path is strictly unbranched. If a branch is encountered, returns None.
        """
        # Consider only carbon neighbors (atomic number 6) excluding the atom we came from.
        neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != coming_from_idx]
        if len(neighbors) == 0:
            return 1  # terminus reached
        elif len(neighbors) == 1:
            next_length = walk_chain(neighbors[0], atom.GetIdx())
            if next_length is None:
                return None
            else:
                return 1 + next_length
        else:
            # Branching encountered; the chain is not a single unbranched fatty acid chain.
            return None
    
    # Start by walking from the unique carbon neighbor of the acid carbon
    chain_from_alpha = walk_chain(carbon_neighbors[0], acid_carbon_idx)
    if chain_from_alpha is None:
        return False, "Fatty acid chain is branched, not a single unbranched linear chain"
    
    # Count total chain length (including the acid carbon)
    total_chain_length = 1 + chain_from_alpha
    
    # Check against threshold: chain length must be >22 to qualify.
    if total_chain_length > 22:
        reason = f"Longest carbon chain from the acid carbon has {total_chain_length} carbons; qualifies as very long-chain fatty acid"
        if total_chain_length > 27:
            reason += " (ultra-long-chain fatty acid)"
        return True, reason
    else:
        return False, f"Longest carbon chain from the acid carbon has {total_chain_length} carbons, which is not >22"

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples:
    test_smiles = "CCCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCC(O)=O"
    result, msg = is_very_long_chain_fatty_acid(test_smiles)
    print(result, msg)