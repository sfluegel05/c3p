"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI entry for 3-hydroxy fatty acids.
Definition: A 3-hydroxy fatty acid (for our purposes) is defined as a molecule that
  (i) contains a free (terminal) carboxylic acid group (COOH),
  (ii) whose alpha carbon (the unique carbon attached to the acid carbon) is linked to at least one candidate beta carbon,
  (iii) and at least one beta carbon (i.e. a carbon 2 bonds away from the acid) bears an –OH group,
  (iv) and the overall contiguous carbon chain (following carbon–carbon bonds) beginning at the acid is of sufficient length (≥5 carbons).
This implementation tries to avoid false positives by ensuring the acid group is terminal and that the longest carbon chain is long.
It relaxes the strict “linear only” rule so that common long-chain fatty acids, even those with minor branches or small rings (e.g. cyclopropyl groups), are accepted.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3‐hydroxy fatty acid based on its SMILES string.
    
    A valid 3‐hydroxy fatty acid (as defined here) must have:
      1. A free, terminal carboxylic acid group (COOH).
      2. The acid carbon is attached to exactly one carbon (the alpha carbon).
      3. Among the alpha carbon's other carbon neighbors (the beta candidates), 
         at least one must have an –OH group attached.
      4. The contiguous carbon chain (connected via C–C bonds starting at the acid carbon)
         has at least 5 carbons.
         
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
    
    # Find terminal carboxylic acid groups.
    # SMARTS for a carboxylic acid: a carbon with one double-bonded O and one -OH.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (COOH) group found"
    
    # Helper: get the longest continuous carbon chain length from a given starting atom.
    # We only follow C-C bonds. 'visited' is used to avoid infinite loops.
    def longest_chain_length(atom, prev_idx, visited):
        best = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:  # only follow carbon neighbors
                continue
            if nbr.GetIdx() == prev_idx:
                continue
            if nbr.GetIdx() in visited:
                continue
            visited.add(nbr.GetIdx())
            length = 1 + longest_chain_length(nbr, atom.GetIdx(), visited)
            visited.remove(nbr.GetIdx())
            if length > best:
                best = length
        return best

    # Examine each carboxylic acid match.
    for match in acid_matches:
        acid_atom = mol.GetAtomWithIdx(match[0])
        
        # Exclude acid groups that are part of a ring.
        if acid_atom.IsInRing():
            continue
        
        # Terminal acid: acid carbon should be bound to exactly one carbon.
        carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue
        alpha_atom = carbon_neighbors[0]
        if alpha_atom is None:
            continue
        # We allow the alpha carbon to be in a ring or branched,
        # but we then examine each candidate beta.
        candidate_beta_atoms = [nbr for nbr in alpha_atom.GetNeighbors() 
                                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != acid_atom.GetIdx()]
        if not candidate_beta_atoms:
            continue
        
        # For at least one candidate beta atom, check that it carries an -OH.
        beta_valid = False
        for beta_atom in candidate_beta_atoms:
            # Check if beta has an -OH attached.
            # For each neighbor of beta, if it is an oxygen (atomic num 8) and has at least one hydrogen attached.
            has_OH = False
            for nbr in beta_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    has_OH = True
                    break
            if has_OH:
                beta_valid = True
                break  # at least one beta candidate qualifies
        if not beta_valid:
            continue
        
        # Compute the longest contiguous carbon chain starting from the acid carbon.
        # We follow C-C bonds. (This includes the acid carbon itself.)
        chain_length = longest_chain_length(acid_atom, None, {acid_atom.GetIdx()})
        if chain_length < 5:
            continue
        
        # All criteria have been met.
        return True, ("Contains a terminal carboxylic acid with a beta-hydroxy substituent "
                      "and a sufficiently long carbon chain (length %d)" % chain_length)
    
    # If no acid group qualified, then the molecule does not meet criteria.
    return False, "No suitable beta-hydroxy fatty acid motif found"

# Example usage:
if __name__ == "__main__":
    # (3R)-3,17-dihydroxymargaric acid (a true positive)
    test_smiles = "OCCCCCCCCCCCCCC[C@@H](O)CC(O)=O"
    result, reason = is_3_hydroxy_fatty_acid(test_smiles)
    print(result, reason)