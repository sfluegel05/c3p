"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI entry for 3-hydroxy fatty acids.
Definition: A 3-hydroxy fatty acid is defined here as a molecule that
  (i) contains a free (terminal) carboxylic acid group,
  (ii) has a beta carbon (the carbon two bonds away from the acid group) that carries an –OH,
  (iii) and presents a sufficiently long, linear (unbranched), acyclic alkyl chain.
This implementation tries to avoid false positives by checking that the acid-group
is truly terminal and that the chain is linear (i.e. no extra carbon branching).
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    
    A valid 3-hydroxy fatty acid (as defined for our purposes) must have:
      1. A free, terminal carboxylic acid group (COOH).
      2. The acid carbon is attached to exactly one carbon (the alpha carbon).
      3. The alpha carbon (adjacent to the acid) is connected to exactly one carbon 
         besides the acid carbon (its beta carbon).
      4. The beta carbon carries an –OH group.
      5. The alkyl chain (starting with the acid and proceeding linearly through alpha,
         beta, and onwards) contains at least 5 carbon atoms.
      6. The carbons in the chain are not parts of rings (i.e. the chain is acyclic).
      7. The chain is linear (i.e. no extra carbon substituents at the key positions).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a 3‐hydroxy fatty acid, 
              False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for a terminal carboxylic acid group.
    # Use a SMARTS that finds the carboxyl carbon: [CX3](=O)[OX2H]
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (COOH) group found"

    # For each matching acid, verify the rest of the criteria.
    for match in acid_matches:
        # match[0] is the carboxyl carbon.
        acid_atom = mol.GetAtomWithIdx(match[0])
        
        # Exclude acid groups that are in rings.
        if acid_atom.IsInRing():
            continue
        
        # For a terminal acid group, the acid carbon should have exactly one carbon neighbor.
        carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue
        
        # Identify the alpha carbon.
        alpha_atom = carbon_neighbors[0]
        if alpha_atom.IsInRing():
            continue
        
        # For a truly linear chain the alpha carbon should have only two carbon neighbors:
        # one being the acid carbon and one being the beta carbon.
        alpha_carbons = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(alpha_carbons) != 2:
            continue
        
        # Identify the beta carbon (the one that is not the acid carbon).
        if alpha_carbons[0].GetIdx() == acid_atom.GetIdx():
            beta_atom = alpha_carbons[1]
        else:
            beta_atom = alpha_carbons[0]
            
        if beta_atom.IsInRing():
            continue
        
        # Check that the beta carbon has an -OH group attached.
        beta_has_OH = False
        for nbr in beta_atom.GetNeighbors():
            # In RDKit the explicit hydrogens are not always present,
            # so we use GetTotalNumHs() to check for attached hydrogen(s).
            if nbr.GetAtomicNum() == 8 and nbr.GetSymbol() == 'O' and nbr.GetTotalNumHs() > 0:
                beta_has_OH = True
                break
        if not beta_has_OH:
            continue

        # Now follow and count a linear carbon chain starting at the acid.
        # Begin with a list: acid_atom -> alpha_atom -> beta_atom.
        chain_atoms = [acid_atom, alpha_atom, beta_atom]
        
        # From beta_atom, extend forward. At each step, there must be exactly one carbon
        # neighbor (other than the atom we just came from) and it must not be in a ring.
        current_atom = beta_atom
        previous_atom = alpha_atom
        while True:
            next_carbons = [nbr for nbr in current_atom.GetNeighbors() 
                            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != previous_atom.GetIdx()]
            if len(next_carbons) == 0:
                # End of chain reached.
                break
            if len(next_carbons) > 1:
                # Branching detected, so the chain is not linear.
                break
            next_atom = next_carbons[0]
            if next_atom.IsInRing():
                break
            chain_atoms.append(next_atom)
            previous_atom = current_atom
            current_atom = next_atom
        
        # Require that the overall chain (from acid to the chain end) has at least 5 carbons.
        if len(chain_atoms) < 5:
            continue
        
        # All criteria are met: classify as a valid 3-hydroxy fatty acid.
        return True, ("Contains a terminal carboxylic acid with a beta-hydroxy substituent "
                      "and a sufficiently long, linear, acyclic carbon chain")
    
    # If no acid group led to the valid pattern, then classify as not a 3-hydroxy fatty acid.
    return False, "No suitable beta-hydroxy fatty acid motif found"

# Example usage (uncomment the following lines to test):
# test_smiles = "OCCCCCCCCCCCCCC[C@@H](O)CC(O)=O"  # (3R)-3,17-dihydroxymargaric acid, a true positive
# result, reason = is_3_hydroxy_fatty_acid(test_smiles)
# print(result, reason)