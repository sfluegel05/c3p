"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: The conjugate base of a fatty acid anion.
Definition: "The conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group of the corresponding fatty acid."
This version refines the search by ensuring:
  1. A deprotonated carboxylate group ([CX3](=O)[O-]) is found.
  2. That carboxylate carbon is terminal (attached to exactly one carbon).
  3. The alpha carbon (the unique neighbour) isn’t overly decorated (its other bonds are only to C and/or O, and total non-carboxylate neighbors ≤2).
  4. The connected carbon chain starting at the carboxylate is sufficiently long (at least MIN_CHAIN_CARBONS) and makes up a good fraction of all carbons.
  
Note: This is a heuristic; one might tweak thresholds (for chain length, fraction, etc.) or further refine the “decoration” rule.
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    Uses several heuristics:
      - It must contain a deprotonated carboxylate group: [CX3](=O)[O-].
      - The carboxylate carbon must be terminal (exactly one carbon neighbour).
      - The unique alpha carbon (carbon neighbour) should have few substituents (other than hydrogen).
      - The contiguous carbon chain (from the carboxylate) must be sufficiently long,
        and represent a large portion of the molecule's carbon framework.
        
    Args:
        smiles (str): SMILES string of the molecule.
    Returns:
        bool: True if the molecule is identified as a fatty acid anion, False otherwise.
        str: A description of the reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # SMARTS for a deprotonated carboxylate group.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No deprotonated carboxylate group found."

    # Helper: find the length of the longest contiguous carbon chain starting at a given carbon atom,
    # following only carbon–carbon bonds.
    def longest_carbon_chain(start_idx, prev_idx, visited):
        current_atom = mol.GetAtomWithIdx(start_idx)
        best_length = 1  # count self
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx == prev_idx or nbr_idx in visited:
                continue
            candidate = 1 + longest_carbon_chain(nbr_idx, start_idx, visited | {nbr_idx})
            if candidate > best_length:
                best_length = candidate
        return best_length

    # Count all carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Settings for our heuristic.
    MIN_CHAIN_CARBONS = 6  # minimum number of carbons expected in a fatty acid chain (including carboxylate carbon)
    MIN_CHAIN_FRACTION = 0.50  # the chain should represent at least 50% of the molecule's carbon atoms

    # Iterate over each match for the carboxylate group.
    for match in mol.GetSubstructMatches(carboxylate_pattern):
        carboxylC_idx = match[0]  # the first atom in SMARTS is the carboxyl carbon.
        carboxylC = mol.GetAtomWithIdx(carboxylC_idx)
        
        # The carboxylate carbon must have exactly one carbon neighbour.
        carbon_neighbors = [nbr for nbr in carboxylC.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue  # not terminal; try next match
        alpha_atom = carbon_neighbors[0]
        alpha_idx = alpha_atom.GetIdx()
        
        # Check the decoration of the alpha carbon:
        # Excluding the carboxylate carbon, collect its other neighbours.
        other_nbrs = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetIdx() != carboxylC_idx]
        # We permit only a limited number of extra substituents (allowing e.g. one hydroxy substituent).
        if len(other_nbrs) > 2:
            continue  # too many substituents; not a simple acyl chain
        
        # Optionally, check that any substituent on the alpha carbon is of an allowed type (C or O).
        for nbr in other_nbrs:
            if nbr.GetAtomicNum() not in (6, 8):
                break
        else:
            # Passed the alpha carbon check.
            # Compute the longest contiguous carbon chain starting from the carboxylate carbon.
            chain_length = longest_carbon_chain(carboxylC_idx, None, {carboxylC_idx})
            if chain_length < MIN_CHAIN_CARBONS:
                return False, f"Alkyl chain too short (only {chain_length} carbons in chain)."
            # Also, ensure that the chain represents a large fraction of the carbon atoms.
            if chain_length < total_carbons * MIN_CHAIN_FRACTION:
                return False, ("While a terminal deprotonated carboxylate exists, " +
                               f"the contiguous carbon chain (length {chain_length}) is not dominant " +
                               f"compared to total carbons ({total_carbons}).")
            return True, "Molecule contains a terminal deprotonated carboxylate group with a sufficiently long, dominant acyl chain."
    
    # If no carboxylate match passes the terminal and chain tests.
    return False, "Carboxylate group found but not in a terminal position typical of a fatty acid anion."
    
# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCCCCC([O-])=O"  # Example: long-chain fatty acid anion (e.g. cerotate)
# print(is_fatty_acid_anion(test_smiles))