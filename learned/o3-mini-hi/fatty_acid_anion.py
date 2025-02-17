"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: The conjugate base of a fatty acid (fatty acid anion)
Definition: "The conjugate base of a fatty acid, arising from deprotonation of the 
carboxylic acid group of the corresponding fatty acid."
Heuristics used in this version:
  1. The molecule must contain at least one deprotonated carboxylate group: [CX3](=O)[O-].
  2. The carboxylate carbon must be terminal – it is attached (via a C–C bond) to exactly one carbon.
  3. From that single (alpha) carbon a contiguous acyl chain is “traced” by following C–C bonds;
     only carbons that are not aromatic are included.
  4. For large molecules (total carbons ≥ 8) the chain must be sufficiently long (at least 6 carbons);
     for small molecules a lower bar is applied.
  5. If more than one terminal carboxylate candidate passes these conditions (for example, in di‐ or poly‐carboxylic species)
     the molecule is not classified as a fatty acid anion.
Note: This heuristic is not perfect. It was tuned so that many accepted fatty acid anions (such as (R)-2-hydroxyhexadecanoate,
prostaglandin H1(1-), etc.) are flagged while many “false positive” cases (where the –COO– function is embedded in polycarboxylates
or aromatic systems) are weeded out.
"""

from rdkit import Chem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is defined as the deprotonated form of a fatty acid – that is,
    it contains a terminal carboxylate group (–C(=O)[O-]) where the carboxylate carbon
    is attached to exactly one carbon (the alpha carbon) from which a single contiguous
    non‐aromatic (aliphatic) chain extends.
    
    For large molecules (total C ≥ 8) the chain must be long (at least 6 C’s). For small molecules,
    this restriction is relaxed to allow examples like 2-hydroxyisobutyrate.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): True and a message if the molecule is classified as a fatty acid anion;
                     False and a message otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Count total carbons in the molecule.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Decide minimum chain length cutoff for a candidate acyl chain.
    # In larger molecules we expect a long-chain fatty acid (≥6 carbons).
    # In smaller molecules, relax the cutoff.
    min_chain = 6 if total_carbons >= 8 else 1
    # Also require that the fatty chain is a dominant feature if possible.
    # (Here we later check the ratio of chain carbons to total carbons.)
    min_chain_ratio = 0.2  # at least 20% of all carbons come from the chain
    
    # Define a SMARTS for a deprotonated carboxylate.
    carboxylate_smarts = "[CX3](=O)[O-]"
    carboxylate_pattern = Chem.MolFromSmarts(carboxylate_smarts)
    matches = mol.GetSubstructMatches(carboxylate_pattern)
    if not matches:
        return False, "No deprotonated carboxylate group found."

    # Helper function: compute the longest contiguous chain (number of carbons)
    # starting from the given atom (by index) along C–C bonds.
    # We restrict to atoms with atomic number 6 that are not in an aromatic ring.
    def longest_chain(atom_idx, prev_idx, visited):
        current_atom = mol.GetAtomWithIdx(atom_idx)
        max_length = 1  # count the current atom
        for nbr in current_atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx == prev_idx or nbr_idx in visited:
                continue
            if nbr.GetAtomicNum() != 6:
                continue
            # Only follow non‐aromatic carbons to avoid chains that lead into rings (like benzene)
            if nbr.GetIsAromatic():
                continue
            # Recurse
            candidate = 1 + longest_chain(nbr_idx, atom_idx, visited | {nbr_idx})
            if candidate > max_length:
                max_length = candidate
        return max_length

    candidate_count = 0
    candidate_reason = ""
    
    # Evaluate each deprotonated carboxylate match.
    for match in matches:
        # In our SMARTS the first atom is the carboxylate carbon.
        carboxyl_c = mol.GetAtomWithIdx(match[0])
        # Find carbon neighbors (exclude oxygens).
        carbon_neighbors = [nbr for nbr in carboxyl_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        # For a terminal group, the carboxylate carbon should be attached to exactly one carbon.
        if len(carbon_neighbors) != 1:
            continue
        alpha_atom = carbon_neighbors[0]
        # Typically in fatty acids the alpha carbon is not inside a ring.
        if alpha_atom.IsInRing():
            continue
        
        # Optionally, one might check that the substituents on alpha are “simple”.
        # (Here we allow any carbon/oxygen attachments and do not further limit the branching.)
        
        # Determine the length of the contiguous aliphatic chain.
        # We start from the alpha carbon and do a DFS along carbon atoms.
        chain_length = longest_chain(alpha_atom.GetIdx(), carboxyl_c.GetIdx(), {alpha_atom.GetIdx()})
        
        # Calculate the ratio of the chain carbons to all carbons.
        chain_ratio = chain_length / total_carbons if total_carbons > 0 else 0
        
        # Check whether the chain is long enough.
        if chain_length >= min_chain or chain_ratio >= min_chain_ratio:
            candidate_count += 1
            candidate_reason = ("Contains a terminal deprotonated carboxylate with an acyl chain of "
                                f"{chain_length} carbons out of {total_carbons} (ratio {chain_ratio:.2f}).")
    
    # If exactly one candidate passes, we classify the molecule as a fatty acid anion.
    if candidate_count == 1:
        return True, candidate_reason
    elif candidate_count == 0:
        return False, "Carboxylate group(s) found but none are attached to a sufficiently long terminal acyl chain."
    else:
        return False, "Multiple terminal carboxylate groups were detected; molecule is not a typical fatty acid anion."

# Example usage:
if __name__ == "__main__":
    # Test with one accepted fatty acid anion:
    test_smiles = "CCCCCCCCCCCCC[C@@H](O)C([O-])=O"  # (R)-2-hydroxyhexadecanoate
    result, reason = is_fatty_acid_anion(test_smiles)
    print(result, reason)