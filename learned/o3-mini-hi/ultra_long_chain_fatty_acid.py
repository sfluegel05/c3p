"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
#!/usr/bin/env python
"""
Classifies: ultra-long-chain fatty acid
Definition: A very long-chain fatty acid is defined as one where the fatty acid portion—
that is, the carbon chain attached to the carboxylic acid (-COOH) group—
consists largely of a single, continuous (acyclic and unbranched) chain whose total number of carbons 
(including the carboxyl carbon) is greater than 27.
"""

from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid.
    
    Strategy:
      1. Parse the SMILES string.
      2. Identify a carboxylic acid group using a SMARTS pattern.
      3. From the carboxyl carbon, find the alpha carbon (a bonded carbon that is not the OH).
      4. Traverse linearly (without branching) from the alpha carbon.
         If at any carbon more than one onward (acyclic) carbon is found (other than the one we came from), 
         we abort because the chain is branched.
      5. Count the number of carbons in the continuous chain (include the carboxyl carbon).
      6. Verify that the chain length > 27.
      7. Also check that nearly all carbons in the molecule belong to the chain.
    
    Args:
      smiles (str): SMILES string of the molecule
      
    Returns:
      bool: True if the molecule qualifies as an ultra-long-chain fatty acid, False otherwise.
      str: Explanation for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find the carboxylic acid group using SMARTS.
    # Pattern matches a trivalent carbon double-bonded to O and single bonded to an OH.
    carboxyl_smarts = "[CX3](=O)[OX2H1]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    matches = mol.GetSubstructMatches(carboxyl_pattern)
    if not matches:
        return False, "Molecule does not contain a carboxylic acid group"
    
    # Assume the first match; the first atom in the match is taken as the carboxyl (carbonyl) carbon.
    carboxyl_idx = matches[0][0]
    carboxyl_atom = mol.GetAtomWithIdx(carboxyl_idx)
    
    # Identify the alpha carbon (the carboxyl carbon most likely has two neighbors:
    # one being the hydroxyl oxygen and the other a carbon from the chain).
    alpha_idx = None
    for nbr in carboxyl_atom.GetNeighbors():
        if nbr.GetAtomicNum() == 6:  # only consider carbon atoms
            alpha_idx = nbr.GetIdx()
            break
    if alpha_idx is None:
        return False, "Carboxyl carbon is not attached to any carbon (no alkyl chain found)"
    
    # Define a function to traverse the chain linearly.
    # starting from a current atom (which is in the chain) with its previous atom,
    # we follow the unique next carbon if it exists.
    def traverse_chain(curr_idx, prev_idx):
        chain_length = 1  # count current carbon
        current_atom = mol.GetAtomWithIdx(curr_idx)
        # Find all carbon neighbors (excluding the previous atom) that are not in rings.
        next_carbons = []
        for nbr in current_atom.GetNeighbors():
            if nbr.GetIdx() == prev_idx:
                continue
            if nbr.GetAtomicNum() == 6 and not nbr.IsInRing():
                next_carbons.append(nbr.GetIdx())
        # If more than one carbon neighbor exists, we assume branching.
        if len(next_carbons) > 1:
            return chain_length, False, "Found branching in the chain"
        # If exactly one exists, continue the linear chain traversal.
        if len(next_carbons) == 1:
            next_idx = next_carbons[0]
            sub_length, linear, reason = traverse_chain(next_idx, curr_idx)
            # Propagate failure of linear chain.
            if not linear:
                return chain_length + sub_length, False, reason
            return chain_length + sub_length, True, ""
        # If none, then termination.
        return chain_length, True, ""
    
    # Traverse starting from the alpha carbon.
    chain_length_from_alpha, is_linear, msg = traverse_chain(alpha_idx, carboxyl_idx)
    # Total chain length includes the carboxyl carbon.
    total_chain_length = 1 + chain_length_from_alpha
    
    if not is_linear:
        return False, f"Chain is branched. {msg}"
    
    # Check if chain is long enough: strictly greater than 27 carbons.
    if total_chain_length <= 27:
        return False, f"Chain length is {total_chain_length} carbons, which is not greater than C27"
    
    # Count the total number of carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # For a true fatty acid, nearly all carbons should be in the acyclic chain.
    # Allow a difference of up to 2 carbons.
    if total_carbons - total_chain_length > 2:
        return False, (f"Chain length is {total_chain_length} carbons but the molecule has {total_carbons} carbons. "
                       "Extra carbon fragments suggest it is not a simple fatty acid.")
    
    return True, f"Chain length is {total_chain_length} carbons, qualifies as ultra-long-chain fatty acid"

# Example usage (uncomment to test):
# if __name__=="__main__":
#     examples = [
#         "CCCCCC\\C=C/CCCCCCCCCC(O)C(O)=O",  # sample ultra-long-chain fatty acid (should be resampled)
#         "OC(CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC(O)=O"  # add other examples if desired
#     ]
#     for smi in examples:
#         res, reason = is_ultra_long_chain_fatty_acid(smi)
#         print(smi, res, reason)