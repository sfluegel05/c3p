"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
#!/usr/bin/env python3
"""
Classifies: long-chain fatty acid (chain length from C13 to C22)
A long-chain fatty acid is defined as a molecule that contains a carboxylic acid group 
whose carbon (the acid carbon) is attached to an unbranched, linear carbon chain 
with a total number of carbons between 13 and 22 (inclusive). 
Note: This is a heuristic method that may fail for molecules with complex structures.
"""

from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A valid long-chain fatty acid must have a carboxylic acid group and an unbranched 
    carbon chain (starting at the acid carbon) whose length (in number of carbons) is between 13 and 22 inclusive.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify carboxylic acid groups. Use a SMARTS pattern matching [CX3](=O)[OX2H]
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_group = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_group)
    
    if not acid_matches:
        return False, "No carboxylic acid group found"
    
    # Try each acid group match to see if it leads to a valid linear chain.
    reasons = []
    for match in acid_matches:
        acid_idx = match[0]  # acid carbon
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        
        # Get all carbon neighbors of the acid carbon.
        acid_carbon_neighbors = [nbr.GetIdx() for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(acid_carbon_neighbors) != 1:
            # If acid carbon is attached to 0 or >1 carbons then it is not connected to a clean, linear alkyl chain.
            reasons.append("Acid carbon does not have exactly one carbon neighbor (found {})".format(len(acid_carbon_neighbors)))
            continue

        # Start the chain from the unique neighbor.
        chain_length = 2  # count acid carbon plus the first chain carbon.
        prev_idx = acid_idx
        current_idx = acid_carbon_neighbors[0]
        branch_detected = False

        # Walk the chain as long as there is exactly one candidate next carbon.
        while True:
            curr_atom = mol.GetAtomWithIdx(current_idx)
            # In our carbon subgraph, list neighbors that are carbons.
            neighbor_idxs = [nbr.GetIdx() for nbr in curr_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
            # Exclude the atom we just came from.
            next_candidates = [idx for idx in neighbor_idxs if idx != prev_idx]
            
            if len(next_candidates) == 0:
                # Terminal carbon reached.
                break
            elif len(next_candidates) > 1:
                # Branching detected: not a strictly linear chain.
                branch_detected = True
                break
            else:
                # Exactly one next carbon: extend the chain.
                chain_length += 1
                prev_idx = current_idx
                current_idx = next_candidates[0]
        
        if branch_detected:
            reasons.append("Branching detected along the chain starting at acid carbon")
            continue
        
        # Now, check if chain_length is within the accepted range.
        if chain_length < 13:
            reasons.append("Longest carbon chain is {} carbons, which is too short".format(chain_length))
            continue
        if chain_length > 22:
            reasons.append("Longest carbon chain is {} carbons, which is too long".format(chain_length))
            continue
        
        # If we reach here then we have a valid unbranched chain of acceptable length.
        return True, "Contains a carboxylic acid group with a connected unbranched carbon chain of {} carbons".format(chain_length)
    
    # If none of the acid groups produced a valid chain, return a failure.
    if reasons:
        return False, reasons[0]
    return False, "No suitable carboxylic acid group with an unbranched carbon chain found"

# Example usage (uncomment to test):
# test_smiles = "OC(=O)CCCCCCCCCCCCCCC([2H])([2H])[2H]"  # Palmitic acid(d3): chain length 16 carbons.
# result, reason = is_long_chain_fatty_acid(test_smiles)
# print(result, reason)