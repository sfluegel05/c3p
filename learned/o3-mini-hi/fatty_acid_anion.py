"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: The conjugate base of a fatty acid (fatty acid anion)
Definition: "The conjugate base of a fatty acid, arising from deprotonation of the 
carboxylic acid group of the corresponding fatty acid."
Heuristics used in this version:
  1. The molecule must contain at least one deprotonated carboxylate group: [CX3](=O)[O-].
  2. The carboxylate carbon must be terminal – it is attached (via a C–C bond) to exactly one carbon (the α–carbon).
  3. From that α–carbon a contiguous chain is “traced” by following only non‐aromatic carbons.
  4. Both an absolute chain‐length and the ratio of chain carbons/total carbons are used; 
     these cutoff values are stricter for larger molecules.
  5. In addition, if the α–carbon has additional substituents (besides the carboxylate and the chain)
     that are not “simple” (only carbon, hydrogen or –OH), the candidate is rejected.
  6. If more than one candidate passes the checks the molecule is not classified as a fatty acid anion.
Note: This heuristic is not perfect but has been tuned based on many true and false positives.
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is defined as the deprotonated form of a fatty acid – that is,
    it contains a terminal carboxylate group (–C(=O)[O-]) in which the carboxyl carbon
    is attached to exactly one carbon (the α–carbon) from which a contiguous aliphatic chain extends.
    
    We compute the longest chain solely along non‐aromatic carbon–carbon bonds and then require that
    (a) the chain is sufficiently long, and (b) the chain carbons represent a high fraction of the overall 
    carbon count. These thresholds depend on the total size of the molecule; very small acids (e.g. 
    2-hydroxyisobutyrate) are treated more leniently.
    
    Additionally, if the α–carbon (the one linked to the carboxylate carbon) has any extra substituents
    other than the carboxylate or carbons (or simple –OH groups), the candidate is rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True with a message if the molecule is classified as a fatty acid anion;
                     False with a message otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Define SMARTS to look for a deprotonated carboxylate group.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not matches:
        return False, "No deprotonated carboxylate group found."
    
    # Helper to get the longest linear chain (list of atom indices) from a starting atom 
    # using DFS. We only follow non-aromatic carbon atoms.
    def longest_chain_path(start_idx, prev_idx, visited):
        current_atom = mol.GetAtomWithIdx(start_idx)
        best_path = [start_idx]
        for nbr in current_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Do not go back to previous or revisit an atom.
            if nbr_idx == prev_idx or nbr_idx in visited:
                continue
            # Only follow carbon atoms.
            if nbr.GetAtomicNum() != 6:
                continue
            # Only follow non‐aromatic carbons so that we do not “jump” into rings.
            if nbr.GetIsAromatic():
                continue
            # Recurse and choose the longest chain 
            new_path = longest_chain_path(nbr_idx, start_idx, visited | {nbr_idx})
            if len(new_path) + 1 > len(best_path):
                best_path = [start_idx] + new_path
        return best_path

    # Allowed elements for substituents on the α–carbon (besides the chain and carboxylate).
    allowed_atomic_nums = {1, 6, 8}  # hydrogen, carbon, oxygen (hydroxyls allowed)
    
    candidate_count = 0
    candidate_reason = ""
    
    # Process each match that is a candidate carboxylate.
    for match in matches:
        # Our SMARTS "[CX3](=O)[O-]" returns (0,1,2) but we consider only the carboxyl carbon:
        carboxyl_c = mol.GetAtomWithIdx(match[0])
        # Get its neighboring carbons.
        carbon_neighbors = [nbr for nbr in carboxyl_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # If not exactly one carbon neighbor the group is not terminal.
        if len(carbon_neighbors) != 1:
            continue
        alpha_atom = carbon_neighbors[0]
        
        # Reject if the α–carbon is in a ring (typical fatty acids develop on an open chain).
        if alpha_atom.IsInRing():
            continue
        
        # Compute the longest contiguous chain starting at the α–carbon.
        chain_path = longest_chain_path(alpha_atom.GetIdx(), carboxyl_c.GetIdx(), {alpha_atom.GetIdx()})
        chain_length = len(chain_path)
        chain_ratio = chain_length / total_carbons if total_carbons > 0 else 0
        
        # Now, check the substituents on the α–carbon.
        # The allowed neighbors are the carboxylate carbon and the one that is included in the chain_path.
        # (If there are extra neighbors with heteroatoms (other than H, C, O), we reject this candidate.)
        extra_substituent_found = False
        for nbr in alpha_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            # Skip if this neighbor is the carboxylate.
            if nbr.GetIdx() == carboxyl_c.GetIdx():
                continue
            # Also skip if this neighbor is the one leading into the longest chain path.
            if nbr_idx in chain_path:
                continue
            # If any other neighbor is present and its atomic number is not allowed, reject.
            if nbr.GetAtomicNum() not in allowed_atomic_nums:
                extra_substituent_found = True
                break
        if extra_substituent_found:
            continue

        # Set thresholds that depend on the size of the molecule.
        if total_carbons < 5:
            min_chain = 2
            min_ratio = 0.5
        elif total_carbons <= 10:
            min_chain = 3
            min_ratio = 0.70
        else:
            min_chain = 6
            min_ratio = 0.75

        if chain_length >= min_chain and chain_ratio >= min_ratio:
            candidate_count += 1
            candidate_reason = (f"Contains a terminal deprotonated carboxylate with an acyl chain "
                                f"of {chain_length} carbons out of {total_carbons} (ratio {chain_ratio:.2f}).")

    if candidate_count == 1:
        return True, candidate_reason
    elif candidate_count == 0:
        return False, ("Carboxylate group(s) found but none are attached to a sufficiently long and 'clean' terminal acyl chain.")
    else:
        return False, ("Multiple terminal carboxylate groups were detected; molecule is not a typical fatty acid anion.")

# Example usage:
if __name__ == "__main__":
    # Test one accepted fatty acid anion:
    test_smiles = "CCCCCCCCCCCCC[C@@H](O)C([O-])=O"  # (R)-2-hydroxyhexadecanoate
    result, reason = is_fatty_acid_anion(test_smiles)
    print(result, reason)