"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: The conjugate base of a fatty acid anion.
Definition: "The conjugate base of a fatty acid, arising from deprotonation of the 
carboxylic acid group of the corresponding fatty acid."
This version refines the search using several heuristics:
  1. It must contain a deprotonated carboxylate group: [CX3](=O)[O-].
  2. That carboxylate carbon must be terminal – i.e. it is attached to exactly one 
     carbon atom, and that carbon (the alpha carbon) is not in a ring.
  3. The alpha carbon’s additional substituents must be limited to carbon or oxygen.
  4. The contiguous carbon chain starting at the carboxylate must be sufficiently long,
     with a minimum length (set here to 6 atoms) or, if shorter, must represent a major
     fraction of the molecule’s carbon framework.
Note: These cutoffs (e.g. minimum chain length, chain fraction) are heuristic.
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    The molecule must have a deprotonated (-[O-]) carboxylate group,
    and that carboxylate must be terminal – attached (via its carbon) to exactly one
    carbon that is not in a ring. Furthermore, the contiguous carbon chain starting
    at that carbon should be sufficiently long (at least MIN_CHAIN_CARBONS), or if it isn’t 
    a dominant portion of the molecule then it is not accepted.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is identified as a fatty acid anion, False otherwise.
        str: A textual description of the reasoning.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Define a SMARTS for a deprotonated carboxylate group.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No deprotonated carboxylate group found."
    
    # Count total carbon atoms in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Settings for our heuristic.
    MIN_CHAIN_CARBONS = 6    # require the acyl chain to be at least 6 carbons (typical for a fatty acid)
    MIN_CHAIN_FRACTION = 0.3 # if the chain is not long, it must still be a large fraction of all carbons

    # Helper function to compute the length of the contiguous carbon chain
    # connected to a given start atom following only C–C bonds.
    def longest_carbon_chain(start_idx, prev_idx, visited):
        best = 1  # count the starting atom
        for nbr in mol.GetAtomWithIdx(start_idx).GetNeighbors():
            # Only proceed if neighbor is carbon
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx == prev_idx or nbr_idx in visited:
                continue
            candidate = 1 + longest_carbon_chain(nbr_idx, start_idx, visited | {nbr_idx})
            if candidate > best:
                best = candidate
        return best

    # Iterate over every match for the carboxylate pattern.
    for match in mol.GetSubstructMatches(carboxylate_pattern):
        # Assume the first atom in the SMARTS is the carboxylate carbon.
        carboxyl_c = mol.GetAtomWithIdx(match[0])
        
        # For a terminal carboxylate the carboxyl carbon must have exactly one carbon neighbour.
        carbon_neighbors = [nbr for nbr in carboxyl_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue  # not terminal, try next match
        
        alpha_atom = carbon_neighbors[0]
        # Exclude cases where the alpha carbon is in a ring 
        # (often seen with sugars or porphyrins, which aren’t in our fatty acid class).
        if alpha_atom.IsInRing():
            continue
        
        # Check that the alpha carbon’s other substituents are not highly decorated.
        other_nbrs = [nbr for nbr in alpha_atom.GetNeighbors() if nbr.GetIdx() != carboxyl_c.GetIdx()]
        # Allow only carbons and oxygens.
        if any(nbr.GetAtomicNum() not in (6, 8) for nbr in other_nbrs):
            continue
        # (Optional) One could also require that the number of these extra substituents 
        # is not too many – here we do not enforce a hard maximum.
        
        # Determine the length of the contiguous carbon chain beginning at the carboxylate.
        chain_length = longest_carbon_chain(carboxyl_c.GetIdx(), None, {carboxyl_c.GetIdx()})
        
        # Enforce a minimum chain length.
        if chain_length < MIN_CHAIN_CARBONS:
            # (e.g. very short acid like 2-hydroxybutyrate would be rejected)
            continue
        
        # If the chain represents only a small fraction of the overall carbons,
        # then the acyl chain is not dominant. However, as a special case if the chain
        # is at least MIN_CHAIN_CARBONS (i.e. 6) we accept molecules even if additional 
        # non-chain groups are present.
        if chain_length < total_carbons * MIN_CHAIN_FRACTION and chain_length < 6:
            continue
        
        return True, "Contains a terminal deprotonated carboxylate with a suitably long acyl chain."
        
    return False, "Carboxylate group found but not in a terminal position typical of a fatty acid anion."

# Example usage:
# test_smiles = "CCCCCCCCCCCCC[C@@H](O)C([O-])=O"  # (R)-2-hydroxyhexadecanoate example
# print(is_fatty_acid_anion(test_smiles))