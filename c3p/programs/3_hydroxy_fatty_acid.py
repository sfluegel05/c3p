"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
"""
Classifies: CHEBI entry for 3-hydroxy fatty acids.
Definition: For our purposes a 3-hydroxy fatty acid must have all of the following:
  (i) A free (terminal) carboxylic acid group (COOH) present as a non‐ring moiety.
  (ii) The carboxyl (acid) carbon is bound to exactly one carbon atom (the alpha carbon).
  (iii) At least one of the other carbon neighbors of the alpha carbon (i.e. a candidate beta carbon)
       carries an –OH group. In addition the alpha and beta carbons must be aliphatic (i.e. not aromatic).
  (iv) The contiguous (C–C) carbon chain starting from the acid carbon is sufficiently long (>=5 carbons)
       when following only bonds between non‐aromatic carbons.
This implementation is tuned to avoid classifying non‐fatty acid compounds (even if they do contain an –OH on a beta candidate)
by eliminating candidates with aromatic or highly decorated substituents.
"""

from rdkit import Chem

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3‐hydroxy fatty acid based on its SMILES string.
    
    Criteria:
      1. Contains a terminal carboxylic acid group (SMARTS: [CX3](=O)[OX2H]) where the acid carbon 
         is not in a ring.
      2. The acid carbon is attached to exactly one carbon (the alpha carbon) which must be non‐aromatic.
      3. Among the alpha carbon's other carbon neighbors (the candidate beta atoms), at least one has an –OH 
         group attached and must also be aliphatic.
      4. The contiguous chain of carbons (following C-C bonds between non‐aromatic carbons) beginning at the acid 
         carbon has a length >= 5.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets all of the criteria for a 3‐hydroxy fatty acid, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for a carboxylic acid group (COOH)
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (COOH) group found"
    
    # Helper: recursively compute the longest contiguous chain of carbons starting from a given atom,
    # following only C-C bonds and only traversing carbons that are not aromatic.
    def longest_chain_length(atom, prev_idx, visited):
        best = 1
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:  # only follow carbon
                continue
            # Only follow aliphatic carbons (aromatic ones are not considered part of a fatty acid chain)
            if nbr.GetIsAromatic():
                continue
            if nbr.GetIdx() == prev_idx:
                continue
            if nbr.GetIdx() in visited:
                continue
            visited.add(nbr.GetIdx())
            chain_len = 1 + longest_chain_length(nbr, atom.GetIdx(), visited)
            visited.remove(nbr.GetIdx())
            if chain_len > best:
                best = chain_len
        return best

    # Iterate over each acid group candidate.
    for match in acid_matches:
        acid_atom = mol.GetAtomWithIdx(match[0])
        
        # Exclude acid groups in rings (must be terminal).
        if acid_atom.IsInRing():
            continue
        
        # The acid carbon must be attached to exactly one carbon (the alpha carbon).
        carbon_neighbors = [nbr for nbr in acid_atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue
        alpha_atom = carbon_neighbors[0]
        # The alpha carbon should be aliphatic.
        if alpha_atom.GetIsAromatic():
            continue
        
        # For candidate beta atoms, look among the other carbon neighbors of the alpha carbon.
        candidate_beta_atoms = [nbr for nbr in alpha_atom.GetNeighbors()
                                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != acid_atom.GetIdx()]
        if not candidate_beta_atoms:
            continue
        
        # At least one candidate beta atom must be aliphatic and bear an -OH.
        beta_valid = False
        for beta_atom in candidate_beta_atoms:
            if beta_atom.GetIsAromatic():
                continue
            # Look for an -OH attached to beta_atom.
            has_OH = False
            for nbr in beta_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetTotalNumHs() > 0:
                    has_OH = True
                    break
            if has_OH:
                beta_valid = True
                break
        if not beta_valid:
            continue
        
        # Compute the longest contiguous carbon chain starting from the acid carbon, following only non-aromatic C-C bonds.
        chain_length = longest_chain_length(acid_atom, None, {acid_atom.GetIdx()})
        if chain_length < 5:
            continue
        
        # If all criteria are met for this acid group candidate then classify as a valid 3-hydroxy fatty acid.
        return True, ("Contains a terminal carboxylic acid with a beta-hydroxy substituent and a sufficiently "
                      "long aliphatic carbon chain (length %d)" % chain_length)
    
    # If no acid group candidate meets all requirements.
    return False, "No suitable 3-hydroxy fatty acid motif found"

# Example usage:
if __name__ == "__main__":
    # Test example: (3R)-3,17-dihydroxymargaric acid 
    test_smiles = "OCCCCCCCCCCCCCC[C@@H](O)CC(O)=O"
    result, explanation = is_3_hydroxy_fatty_acid(test_smiles)
    print(result, explanation)