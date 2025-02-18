"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI entry for 3‐hydroxy fatty acids.
Definition: A 3‐hydroxy fatty acid is any fatty acid that has a free, terminal carboxylic acid
group (COOH) with the acid carbon bound to exactly one carbon. Moreover, one of the other
neighbors (the candidate beta, i.e. 3‐position relative to the COOH) must have a free –OH group.
Also, the (linear) contiguous aliphatic chain (defined by unbranched, non‐aromatic C–C bonds)
starting at the acid group must have a minimum length (>=5 carbons).
This implementation attempts to be more stringent than lesser definitions.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3‐hydroxy fatty acid based on its SMILES string.
    
    Criteria:
      1. Contains a terminal (non‐ring) carboxylic acid group (COOH; SMARTS: [CX3](=O)[OX2H]) 
         where the acid carbon has exactly one carbon neighbor.
      2. The unique carbon neighbor (alpha carbon) must be aliphatic.
      3. Among the other carbon neighbors of the alpha carbon (i.e. candidate beta atoms),
         at least one must be aliphatic and bear an –OH group.
         (This –OH is taken to indicate a beta hydroxy substituent.)
      4. When “walking” from the acid carbon along carbon–carbon bonds in a linear (unbranched)
         fashion (only following non‐aromatic carbons), the chain must have a length >= 5.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if all criteria are met, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for a carboxylic acid group (free COOH)
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (COOH) group found"
    
    # Helper: Walk the chain linearly from a starting carbon.
    # We only follow C–C bonds where the current carbon (except the starting node)
    # has exactly one new carbon neighbor (i.e. no branching).
    def linear_chain_length(start_atom, coming_from=None):
        length = 1  # count the start_atom itself
        current = start_atom
        prev_idx = coming_from.GetIdx() if coming_from is not None else -1
        while True:
            # get neighbors that are carbons, not aromatic and not the one we came from.
            nbrs = [nbr for nbr in current.GetNeighbors() 
                    if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic() 
                    and nbr.GetIdx() != prev_idx]
            if len(nbrs) == 1:
                length += 1
                prev_idx = current.GetIdx()
                current = nbrs[0]
            else:
                # if branching or dead end, stop the linear walk.
                break
        return length

    # Iterate over each candidate COOH group
    for match in acid_matches:
        acid_atom = mol.GetAtomWithIdx(match[0])
        # Exclude if the acid group is inside a ring.
        if acid_atom.IsInRing():
            continue
        
        # The acid carbon should have exactly one carbon neighbor (the alpha carbon).
        carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue
        alpha_atom = carbon_neighbors[0]
        if alpha_atom.GetIsAromatic():
            continue
        
        # Now, among the alpha carbon's neighbors (other than the acid atom), we check for a candidate beta.
        beta_found = False
        for nbr in alpha_atom.GetNeighbors():
            # Only consider carbon neighbors (the beta candidate) and exclude the acid group.
            if nbr.GetAtomicNum() != 6 or nbr.GetIdx() == acid_atom.GetIdx():
                continue
            if nbr.GetIsAromatic():
                continue
            # Check if this beta candidate has an –OH substituent.
            # Look for an oxygen attached to beta that is not doubly bound (i.e. bearing H(s)).
            for beta_nbr in nbr.GetNeighbors():
                if beta_nbr.GetAtomicNum() == 8 and beta_nbr.GetTotalNumHs() > 0:
                    beta_found = True
                    break
            if beta_found:
                break
                    
        if not beta_found:
            continue  # no qualifying beta hydroxy substituent on this acid candidate
        
        # Now check the aliphatic chain length.
        # We start from the acid atom and walk along its only carbon neighbor (alpha)
        chain_length = linear_chain_length(acid_atom, coming_from=None)
        if chain_length < 5:
            continue

        # All criteria met.
        return True, ("Contains a terminal carboxylic acid with a beta-hydroxy substituent and a sufficiently "
                      "long linear aliphatic carbon chain (chain length %d)" % chain_length)

    # If no candidate COOH met the full requirements.
    return False, "No suitable 3-hydroxy fatty acid motif found"

# Example usage (this block can be removed or commented out if imported as a module):
if __name__ == "__main__":
    examples = [
        ("OCCCCCCCCCCCCCC[C@@H](O)CC(O)=O", "(3R)-3,17-dihydroxymargaric acid"),
        ("C(C(O)=O)C(CCCC)O", "3-hydroxyheptanoic acid"),
        ("OCCCCCCCCCCCCCCCCCCCC[C@@H](O)CC(O)=O", "(3R)-3,22-dihydroxytricosanoic acid"),
        ("OC(C[C@@H](CC)O)=O", "(R)-3-hydroxypentanoic acid"),
        ("OCCCCCC[C@@H](O)CC(O)=O", "(3R)-3,9-dihydroxynonanoic acid"),
    ]
    for smi, name in examples:
        valid, reason = is_3_hydroxy_fatty_acid(smi)
        print(f"{name}: {valid}. Reason: {reason}")