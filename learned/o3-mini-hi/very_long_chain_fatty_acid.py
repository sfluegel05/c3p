"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid.
Heuristic:
  1. Molecule must contain a free (protonated) carboxylic acid group (C(=O)O).
  2. The acid carbon (first atom in the SMARTS match) must be attached to exactly one carbon—a chain initiator.
  3. From that neighbor, we “walk” atom‐by‐atom. At each step the current chain carbon (excluding the atom we came from)
     must have exactly one new carbon neighbor. If there are zero such neighbors the chain is terminated,
     while if there is more than one then the chain is branched and the molecule is rejected.
  4. Each chain atom is required not to be in a ring.
  5. Finally, counting the acid carbon the chain must contain >22 carbons (or >27 for ultra‐long‐chain fatty acids).
"""

from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    
    A very long-chain fatty acid is heuristically defined as a molecule that contains a single free
    carboxylic acid group (C(=O)O). In the chosen acid motif the acid carbon (first atom in the SMARTS)
    should be attached to exactly one carbon atom, which starts the main unbranched acyclic chain.
    We then follow the chain atom-by-atom, at each step requiring that the current chain carbon
    (ignoring the previous atom) has exactly one new carbon neighbor (no extra direct C–C substituents).
    The chain ends when no new carbon neighbor can be found. The total number of carbons (including the acid carbon)
    must be greater than 22 (with >27 noted as ultra-long-chain).
    
    Args:
        smiles (str): SMILES string for the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a very long-chain fatty acid, False otherwise.
        str: A reason message for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for a free (protonated) carboxylic acid group: C(=O)O.
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No free (protonated) carboxylic acid group found"
    
    # Use the first matching acid group.
    # In the SMARTS the acid carbon is the first atom.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # For the acid carbon, determine its carbon neighbors.
    acid_c_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(acid_c_neighbors) != 1:
        return False, f"Acid carbon has {len(acid_c_neighbors)} carbon neighbors; expected exactly one chain initiator"
    
    # Initialize the chain. We'll consider the acid carbon as the start.
    chain_length = 1  # acid carbon count
    # Set the first link in the chain.
    prev_atom = acid_carbon
    current_atom = acid_c_neighbors[0]
    chain_length += 1

    # Walk along the chain.
    while True:
        # Ensure the current chain atom is not in any ring.
        if current_atom.IsInRing():
            return False, f"Chain atom with idx {current_atom.GetIdx()} is in a ring; chain must be acyclic"
        
        # Look for carbon neighbors of current_atom excluding the atom we just came from.
        next_candidates = [
            nbr for nbr in current_atom.GetNeighbors() 
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx()
        ]
        
        if len(next_candidates) == 0:
            # Chain terminated.
            break
        if len(next_candidates) > 1:
            return False, f"Fatty acid chain is branched at atom idx {current_atom.GetIdx()}"
        
        # Exactly one new carbon found; extend the chain.
        next_atom = next_candidates[0]
        chain_length += 1
        prev_atom, current_atom = current_atom, next_atom

    if chain_length <= 22:
        return False, f"Longest carbon chain (length {chain_length}) is not >22"
    
    reason = f"Longest carbon chain from the acid carbon has {chain_length} carbons; qualifies as very long-chain fatty acid"
    if chain_length > 27:
        reason += " (ultra-long-chain fatty acid)"
    return True, reason

# Example usage (this block is self-contained):
if __name__ == "__main__":
    # Test example: you may replace this SMILES with any from the provided list.
    test_smiles = "CCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC(O)=O".replace(" ", "")
    result, msg = is_very_long_chain_fatty_acid(test_smiles)
    print(result, msg)