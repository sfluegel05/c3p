"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: very long-chain fatty acid.
Heuristic:
  1. The molecule must contain a free (protonated) carboxylic acid group (SMARTS: C(=O)O).
  2. The acid carbon (first atom in the pattern) must have exactly one carbon neighbor.
  3. From that neighbor, follow a chain atom-by-atom. At each step, the carbon (ignoring the one we came from)
     must be unique (i.e. no branches) and not be in a ring.
  4. In addition, no chain atom may have extra carbon substituents attached directly.
  5. Finally, the total number of carbons in the main chain (including the acid carbon) must be >22
     (long-chain) and >27 is noted as ultra-long-chain.
     
Note: This heuristic does allow a small substituent (if linked via an intervening heteroatom, e.g. a methoxy)
attached to the chain. However, any extra direct (C–C) substituents not incorporated in the main unbranched chain
will cause the molecule to be rejected.
"""

from rdkit import Chem

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    
    A very long-chain fatty acid is heuristically defined as a molecule containing a single free (protonated)
    carboxylic acid group (C(=O)O) whose acid carbon is attached to a single, unbranched acyclic chain.
    This chain is “walked” atom-by-atom by requiring that each chain carbon (after the acid carbon) 
    has exactly one other carbon neighbor (aside from the previous atom). Furthermore, if any chain atom 
    has an extra directly attached carbon (i.e. a branch) the molecule is rejected – note that groups such as 
    methoxy are allowed since the chain atom is attached to oxygen and then to a carbon.
    
    Finally, the number of carbons in the main chain (counting the acid carbon) must be >22. If >27 we note it as 
    ultra-long-chain.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a very long-chain fatty acid, False otherwise.
        str: Reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Only consider the free (protonated) carboxylic acid: C(=O)O
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No free (protonated) carboxylic acid group found"

    # Use first matching acid group.
    # In the SMARTS the acid carbon appears as the first atom.
    acid_carbon_idx = acid_matches[0][0]
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    
    # The acid carbon should have exactly one carbon neighbor (the start of the chain)
    c_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(c_neighbors) != 1:
        return False, f"Acid carbon has {len(c_neighbors)} carbon neighbors; expected exactly one chain initiation"
    
    chain_indices = set()
    chain_indices.add(acid_carbon_idx)
    chain_atoms = []  # list of atoms along the chain (in order)
    chain_atoms.append(acid_carbon)
    
    # Start with the unique chain neighbor.
    current_atom = c_neighbors[0]
    chain_indices.add(current_atom.GetIdx())
    chain_atoms.append(current_atom)
    chain_length = 2  # acid carbon + first neighbor

    prev_idx = acid_carbon_idx
    
    # Walk the chain as long as exactly one carbon neighbor (other than the atom we came from) is present.
    while True:
        # Check that current atom is not in any ring.
        if current_atom.IsInRing():
            return False, f"Chain atom with idx {current_atom.GetIdx()} is in a ring; chain must be acyclic"
        
        # Check for extra direct carbon substituents on the current chain atom.
        # (Allow the atom we came from; all other C neighbors not in the main chain are not allowed.)
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in chain_indices:
                return False, f"Extra carbon substituent detected on chain atom idx {current_atom.GetIdx()}; chain is branched"
        
        # Get carbon neighbors (exclude what we came from)
        nbr_carbons = [nbr for nbr in current_atom.GetNeighbors() 
                       if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_idx]
        
        if len(nbr_carbons) == 0:
            # Chain terminated
            break
        if len(nbr_carbons) > 1:
            return False, "Fatty acid chain is branched, not a single unbranched linear chain"
        
        # Continue the chain.
        next_atom = nbr_carbons[0]
        chain_length += 1
        chain_indices.add(next_atom.GetIdx())
        chain_atoms.append(next_atom)
        
        # Step along the chain.
        prev_idx = current_atom.GetIdx()
        current_atom = next_atom

    # As an additional purity filter, nearly all carbon atoms of the molecule should be in the main chain.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    extra_carbons = total_carbons - len(chain_indices)
    # Allow up to 1 extra carbon (this may cover minor substituents attached via heteroatoms)
    if extra_carbons > 1:
        return False, f"Too many extra carbons outside the main chain ({extra_carbons} extra carbons detected)"
    
    if chain_length <= 22:
        return False, f"Longest carbon chain from the acid carbon has {chain_length} carbons, which is not >22"
    
    reason = f"Longest carbon chain from the acid carbon has {chain_length} carbons; qualifies as very long-chain fatty acid"
    if chain_length > 27:
        reason += " (ultra-long-chain fatty acid)"
    
    return True, reason

# Example usage (this code block is self-contained):
if __name__ == "__main__":
    # A test example from the provided list. You can replace this with other SMILES strings.
    test_smiles = "CCCC C\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC(O)=O".replace(" ", "")
    result, msg = is_very_long_chain_fatty_acid(test_smiles)
    print(result, msg)