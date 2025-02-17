"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid (chain length > C22; > C27 = ultra-long-chain)
Heuristic: The molecule must have a carboxylic acid group. Then, starting from the acid carbon,
we identify the unique carbon neighbor (the alpha carbon) that starts an unbranched (linear) carbon chain.
If the chain (including the acid carbon) has more than 22 carbons in succession, then
the molecule qualifies as a very long-chain fatty acid.
"""
from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as having an unbranched fatty acid chain (starting
    from the carboxylic acid carbon) with more than 22 contiguous carbon atoms.
    (Fatty acids with more than 27 chain carbons are typically called ultra-long-chain.)
    
    The algorithm:
      1. Confirms that a carboxylic acid group (protonated or deprotonated) is present.
      2. Identifies the acid's carbon, then selects its unique carbon neighbor (the alpha-carbon).
         If there is not exactly one such neighbor, we cannot clearly define a fatty acid chain.
      3. “Walks” down the chain: from each current carbon, if exactly one new carbon neighbor 
         (excluding the one we came from) is present, add it into the chain.
      4. Returns True if the chain length (counting the acid carbon as the first in the chain)
         is >22; otherwise returns False.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule qualifies as a very long-chain fatty acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a carboxylic acid group.
    # We consider both protonated and deprotonated forms.
    acid_smarts1 = Chem.MolFromSmarts("C(=O)[OH]")
    acid_smarts2 = Chem.MolFromSmarts("C(=O)[O-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts1) + mol.GetSubstructMatches(acid_smarts2)
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Use the first match: the first atom in the match is the acid carbon.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # Among acid carbon's neighbors, select only carbon atoms.
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, f"Expected exactly one carbon neighbor off the acid carbon, found {len(carbon_neighbors)}"
    
    # Start the chain with the acid carbon and its unique carbon neighbor.
    chain_indices = [acid_carbon_idx]
    current_atom = carbon_neighbors[0]
    chain_indices.append(current_atom.GetIdx())
    previous_atom_idx = acid_carbon_idx

    # Walk down the chain as long as the chain is linear 
    # (i.e. each carbon (except possibly termini) should have exactly one new carbon neighbor).
    while True:
        # For the current carbon atom get its carbon neighbors, excluding the one we came from.
        nbrs = [nbr for nbr in current_atom.GetNeighbors() 
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != previous_atom_idx]
        # If exactly one new carbon is present, then continue along the chain.
        if len(nbrs) == 1:
            next_atom = nbrs[0]
            chain_indices.append(next_atom.GetIdx())
            previous_atom_idx = current_atom.GetIdx()
            current_atom = next_atom
        else:
            # If none or more than one are found, we have reached the chain terminus 
            # (or the chain is branched) so we stop.
            break

    chain_length = len(chain_indices)
    
    # A valid fatty acid should have an unbranched chain.
    if chain_length <= 2:
        return False, f"Chain length of {chain_length} (acid carbon + {chain_length-1} in chain) is too short to be considered a fatty acid"

    # Check chain length against threshold.
    if chain_length <= 22:
        return False, f"Longest carbon chain from the acid carbon has {chain_length} carbons, which is not >22"
    else:
        reason = f"Longest carbon chain from the acid carbon has {chain_length} carbons; qualifies as very long-chain fatty acid"
        if chain_length > 27:
            reason += " (ultra-long-chain fatty acid)"
        return True, reason

# Example usage:
if __name__ == "__main__":
    # Test with one of the examples: 
    test_smiles = "CCCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCC(O)=O"
    result, msg = is_very_long_chain_fatty_acid(test_smiles)
    print(result, msg)