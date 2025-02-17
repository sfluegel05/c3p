"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid (chain length > C22; > C27 = ultra-long-chain)
Heuristic improvements:
  1. Ensure a valid carboxylic acid group is present (protonated or deprotonated).
  2. For the acid carbon, only one carbon neighbor is allowed.
  3. Starting from that neighbor, follow the chain atom–by–atom. At each step, the carbon (ignoring the atom we came from)
     must have exactly one carbon neighbor. (That is, no extra branches.)
  4. Also require that none of the chain carbons lie in a ring.
  5. As an extra filter, make sure that nearly all carbon atoms in the molecule belong to the chain (allowing a small allowance for minor substituents,
     such as a methoxy group at C2).
  6. Finally, count the chain length (including the acid carbon); if >22 then the molecule qualifies as a very long‐chain fatty acid.
     If >27, note it is ultra-long‐chain.
     
Note: In fatty acid nomenclature the “chain length” is often taken to be the total number of carbons (including the acid carbon).
This heuristic is not perfect but improves on the previous one.
"""

from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    
    A very long-chain fatty acid is (heuristically) defined as a molecule having a carboxylic acid group
    (protonated [C(=O)O] or deprotonated [C(=O)[O-]]) where the acid carbon is attached to a single unbranched, 
    acyclic (non-ring) linear chain of carbons. To be classified, the total number of carbons in the main chain
    (starting from the acid carbon) must be greater than 22. If the chain is even longer (>27 carbons), 
    it is noted as ultra-long-chain.
    
    In addition to tracing the chain, an extra “purity” filter is applied:
    nearly all carbon atoms in the molecule (all but up to two) should lie in the main chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a very long-chain fatty acid, False otherwise.
        str: Reason for the decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS for carboxylic acid: allow both protonated and deprotonated forms.
    acid_pattern1 = Chem.MolFromSmarts("C(=O)[OH]")
    acid_pattern2 = Chem.MolFromSmarts("C(=O)[O-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern1) + mol.GetSubstructMatches(acid_pattern2)
    if not acid_matches:
        return False, "No carboxylic acid group found"

    # For now, we use the first found acid group.
    # In the SMARTS match the first atom is the acid C.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)

    # Among acid carbon's neighbors, select only carbon atoms.
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, f"Expected exactly one carbon neighbor off the acid carbon, found {len(carbon_neighbors)}."

    # Now, walk the potential fatty acid chain starting from that unique neighbor.
    # We also record the set of carbon atoms that are part of the main chain.
    chain_indices = {acid_carbon_idx}  # include the acid carbon
    current_atom = carbon_neighbors[0]
    chain_indices.add(current_atom.GetIdx())
    chain_length = 2  # count acid carbon and its first neighbor

    # Walk while the chain is unbranched, and check that no chain atom is in a ring.
    prev_idx = acid_carbon_idx
    while True:
        # Check that current atom is not in any ring.
        if current_atom.IsInRing():
            return False, "Chain atom is in a ring; fatty acid chain must be acyclic."

        # Get current carbon neighbors excluding the one we just came from.
        nbr_carbons = [nbr for nbr in current_atom.GetNeighbors()
                       if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_idx]
        # If no further carbon neighbor, chain terminates.
        if len(nbr_carbons) == 0:
            break
        # If more than one, then a branch is encountered.
        if len(nbr_carbons) > 1:
            return False, "Fatty acid chain is branched, not a single unbranched linear chain."
        # Otherwise, exactly one; continue the chain.
        next_atom = nbr_carbons[0]
        chain_indices.add(next_atom.GetIdx())
        chain_length += 1
        prev_idx = current_atom.GetIdx()
        current_atom = next_atom

    # Now, use the heuristic that for a free fatty acid almost all carbon atoms in the molecule
    # should be in the main chain. Allow a small difference (allow up to two extra carbons)
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    extra_carbons = total_carbons - len(chain_indices)
    if extra_carbons > 2:
        return False, (f"Too many carbon atoms outside the main chain "
                       f"({extra_carbons} extra carbons detected), so unlikely to be a free fatty acid.")

    # Check chain length against the threshold.
    # Here the chain length includes the acid carbon.
    if chain_length > 22:
        reason = f"Longest carbon chain from the acid carbon has {chain_length} carbons; qualifies as very long-chain fatty acid"
        if chain_length > 27:
            reason += " (ultra-long-chain fatty acid)"
        return True, reason
    else:
        return False, f"Longest carbon chain from the acid carbon has {chain_length} carbons, which is not >22."


# Example usage (this code block is self-contained):
if __name__ == "__main__":
    # Test one example from the provided list:
    test_smiles = "CCCCCC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCCCCC(O)=O"
    result, msg = is_very_long_chain_fatty_acid(test_smiles)
    print(result, msg)