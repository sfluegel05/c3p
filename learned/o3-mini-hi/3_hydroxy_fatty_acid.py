"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI entry for 3-hydroxy fatty acids.
Definition: Any fatty acid with a hydroxy functional group in the beta- (3-) position.
This implementation checks that the molecule contains a terminal (free) carboxylic acid 
group, that the carbon directly alpha to the acid is linked to a beta carbon bearing an –OH,
and that the attached alkyl chain is sufficiently long and linear.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A valid 3-hydroxy fatty acid must have:
      1. A free (terminal) carboxylic acid group (COOH).
      2. The carbon directly alpha to the acid must be connected to a beta carbon 
         that carries an –OH group.
      3. The aliphatic chain (from the acid group onward) should be sufficiently long,
         suggesting a fatty acid rather than a fragment of a larger (e.g., polycyclic) structure.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is a 3-hydroxy fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Attempt to parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # First, search for a free carboxylic acid group.
    # The SMARTS below looks for a carbon (sp2) with a double bond to O and a single bond to an OH.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (COOH) group found"
    
    # For each carboxylic acid match, try to verify the beta-hydroxy fatty acid motif.
    for match in acid_matches:
        # In this SMARTS, match[0] is the acid carbon.
        acid_idx = match[0]
        acid_atom = mol.GetAtomWithIdx(acid_idx)
        
        # For a terminal (free) acid group, the acid carbon should be attached to exactly one carbon.
        # (The other bond goes to the -OH.)
        carbon_neighbors = [n for n in acid_atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            # Skip if acid carbon is not terminal (e.g. if it is connected to more than one carbon).
            continue
        
        # Also, if the acid carbon is in a ring, it likely does not belong to a fatty acid.
        if acid_atom.IsInRing():
            continue
        
        # Identify the alpha carbon (the carbon bonded to the acid carbon).
        alpha_atom = carbon_neighbors[0]
        
        # Look for the beta carbon. For a conventional fatty acid, the alpha carbon should have 
        # one carbon neighbor besides the acid carbon.
        beta_candidates = [n for n in alpha_atom.GetNeighbors() if n.GetIdx() != acid_idx and n.GetAtomicNum() == 6]
        if not beta_candidates:
            continue  # no beta carbon found, so move to next acid match.
        
        # Check whether any beta candidate carries a hydroxy group (–OH).
        beta_has_OH = False
        for beta_atom in beta_candidates:
            for nb in beta_atom.GetNeighbors():
                # Look for oxygen atoms that have at least one hydrogen attached.
                if nb.GetAtomicNum() == 8 and nb.GetTotalNumHs() > 0:
                    beta_has_OH = True
                    break
            if beta_has_OH:
                break
        if not beta_has_OH:
            continue  # this acid group does not have a beta-OH arrangement.
        
        # Next, check that the alkyl chain is long enough to be considered a fatty acid.
        # We perform a simple depth-first search (DFS) from the alpha carbon (excluding the acid)
        # along carbon atoms that are not in rings.
        def dfs_chain(atom, prev_idx, count, visited):
            max_len = count
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != prev_idx:
                    # Skip if the carbon is in a ring; we favor flexible acyclic chains.
                    if neighbor.IsInRing():
                        continue
                    if neighbor.GetIdx() in visited:
                        continue
                    visited.add(neighbor.GetIdx())
                    new_len = dfs_chain(neighbor, atom.GetIdx(), count + 1, visited)
                    visited.remove(neighbor.GetIdx())
                    if new_len > max_len:
                        max_len = new_len
            return max_len
        
        visited = set([acid_atom.GetIdx(), alpha_atom.GetIdx()])
        chain_length = dfs_chain(alpha_atom, acid_atom.GetIdx(), 1, visited)
        
        # Require a minimum chain length; here 6 (alpha plus five more carbons) is set arbitrarily.
        if chain_length < 6:
            continue
        
        # If all criteria are met, classify the molecule as a valid 3-hydroxy fatty acid.
        return True, ("Contains a terminal carboxylic acid with a beta-hydroxy substituent "
                      "on a sufficiently long, acyclic aliphatic chain")
    
    # If no acid group led to a valid disposition, then classify as not a 3-hydroxy fatty acid.
    return False, "No suitable beta-hydroxy fatty acid motif found"

# Example (uncomment to test):
# test_smiles = "OCCCCCCCCCCCCCC[C@@H](O)CC(O)=O"  # (3R)-3,17-dihydroxymargaric acid
# result, reason = is_3_hydroxy_fatty_acid(test_smiles)
# print(result, reason)