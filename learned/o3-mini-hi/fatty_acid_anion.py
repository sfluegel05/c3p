"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: The conjugate base of a fatty acid (i.e. fatty acid anion)
Definition: "The conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group of the corresponding fatty acid."
This version refines the search by ensuring:
  1. A deprotonated carboxylate group ([CX3](=O)[O-]) is found.
  2. The carboxylate carbon is terminal (has exactly one carbon neighbour).
  3. The alpha carbon (that unique carbon neighbour) is not heavily decorated with heteroatoms (only C or O allowed besides the carboxylate).
  4. The contiguous carbon chain length (starting from the carboxylate) is at least 4 carbons.
  
Note: This is a heuristic; one might refine the chain rules or thresholds.
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    Uses several heuristics:
      - It must contain a deprotonated carboxylate group: [CX3](=O)[O-].
      - This carboxylate group must be terminal: its carbon must be attached to exactly one carbon (the alpha carbon).
      - The alpha carbon should only be bonded (besides to the carboxylate) to carbons or possibly an oxygen (to allow hydroxy variants).
      - The contiguous carbon chain (starting at the carboxylate and including the alpha carbon and its extension) is required to be at least 4 carbons long.
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if molecule is identified as a fatty acid anion, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # SMARTS for deprotonated carboxylate.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No deprotonated carboxylate group found."
    
    # Helper: recursively find the maximum number of carbon atoms in a simple (acyclic) path
    def longest_carbon_chain(atom_idx, prev_idx, visited):
        # visited: set of atom indices in the current path.
        current_atom = mol.GetAtomWithIdx(atom_idx)
        best = 1  # count self
        for nbr in current_atom.GetNeighbors():
            # Only follow carbon neighbors (atomic number 6)
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx == prev_idx or nbr_idx in visited:
                continue
            candidate = 1 + longest_carbon_chain(nbr_idx, atom_idx, visited | {nbr_idx})
            if candidate > best:
                best = candidate
        return best

    MIN_TOTAL_CARBONS = 4  # minimum carbons in fatty acid (including the carboxylate carbon)
    
    # Iterate over all carboxylate group matches.
    for match in mol.GetSubstructMatches(carboxylate_pattern):
        # According to our SMARTS the first atom is the carboxyl carbon.
        carboxylC_idx = match[0]
        carboxylC = mol.GetAtomWithIdx(carboxylC_idx)
        
        # Determine the carbon neighbors of the carboxyl carbon.
        carbon_neighbors = [nbr for nbr in carboxylC.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            # If not exactly one carbon neighbour, then this carboxylate is not terminal.
            continue
        
        alpha_atom = carbon_neighbors[0]
        alpha_idx = alpha_atom.GetIdx()
        
        # Look at the other neighbors of the alpha carbon (excluding the carboxyl carbon)
        for nbr in alpha_atom.GetNeighbors():
            if nbr.GetIdx() == carboxylC_idx:
                continue
            # Allow carbon and oxygen (to allow hydroxy substitution)
            if nbr.GetAtomicNum() not in (6, 8):
                # This alpha group is too decorated – likely not a simple acyl chain.
                break
        else:
            # Passed the alpha carbon heteroatom check.
            # Compute the length of this continuous carbon chain.
            # We do a DFS starting at the carboxylate carbon, but only follow carbon bonds.
            total_chain = longest_carbon_chain(carboxylC_idx, None, {carboxylC_idx})
            if total_chain < MIN_TOTAL_CARBONS:
                # Chain is too short.
                return False, f"Alkyl chain too short (only {total_chain} carbons in chain)."
            return True, "Molecule contains a terminal deprotonated carboxylate group with a sufficiently long acyl chain."
    
    # If no match satisfies the terminal condition.
    return False, "Carboxylate group found but not in a terminal position typical of a fatty acid anion."

# Example usage (uncomment the following lines to test):
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCC([O-])=O"  # cerotate example (true positive)
# print(is_fatty_acid_anion(test_smiles))