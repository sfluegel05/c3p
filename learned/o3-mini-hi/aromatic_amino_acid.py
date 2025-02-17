"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acid
Definition: An amino acid whose structure includes an aromatic ring. In other words,
the molecule must contain an integrated (free) alpha–amino acid backbone (with a free –NH2 and a free –COOH on the same carbon)
and that backbone must carry an aromatic side chain.
"""

from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    For a positive classification, the molecule must contain an alpha–amino acid backbone
    (an sp3 alpha–carbon attached to a free (non‐amidic) amino group and the directly connected carboxyl carbon)
    and the side chain attached to the alpha–carbon must include an aromatic ring.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an aromatic amino acid, False otherwise.
        str: Explanation of the classification result.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for better matching.
    mol = Chem.AddHs(mol)
    
    # --- Step 1. Look for an alpha–amino acid backbone ---
    # The SMARTS is intended to match the alpha–carbon (non‐aromatic, sp3 with exactly one hydrogen),
    # a free amino group, and the carboxyl group (we match the carboxyl carbon and one carboxyl oxygen).
    backbone_smarts = "[C;!a;H1](N)(C(=O)[O;H1,-])"
    aa_backbone = Chem.MolFromSmarts(backbone_smarts)
    matches = mol.GetSubstructMatches(aa_backbone)
    if not matches:
        return False, "Missing integrated alpha–amino acid backbone (free amino and acid groups on same carbon)"
    
    # --- Step 2. Check for aromatic side chain ---
    # The typical match will produce a tuple with 4 indices:
    #   index 0: alpha–carbon
    #   index 1: free amino nitrogen
    #   index 2: carboxyl carbon (directly attached to alpha–carbon)
    #   index 3: carboxyl oxygen (attached to the carboxyl carbon)
    # We ignore the oxygen (match[3]) and focus on the neighbors of the alpha–carbon.
    
    def side_chain_has_aromatic(start_idx, visited, max_depth=6):
        """
        Depth-first search starting from atom with index start_idx to detect an aromatic atom.
        Traverses bonds up to a specified max_depth.
        """
        if max_depth < 0:
            return False
        atom = mol.GetAtomWithIdx(start_idx)
        if atom.GetIsAromatic():
            return True
        for nbr in atom.GetNeighbors():
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            visited.add(nbr_idx)
            if side_chain_has_aromatic(nbr_idx, visited, max_depth - 1):
                return True
        return False

    aromatic_sidechain_found = False
    # Iterate over each backbone match.
    for match in matches:
        if len(match) < 3:
            continue  # safety check; we expect 3 or 4 members
        alpha_idx = match[0]
        amine_idx = match[1]
        carboxyl_c_idx = match[2]  # the carboxyl carbon directly attached to the alpha–carbon

        # Get all neighbors of the alpha–carbon.
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        neighbor_indices = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors()]
        # Exclude the amino nitrogen and the carboxyl carbon.
        side_chain_indices = [idx for idx in neighbor_indices if idx not in (amine_idx, carboxyl_c_idx)]
        
        # If there is no side chain (e.g., glycine), continue checking other matches.
        if not side_chain_indices:
            continue
        
        # For each candidate side-chain atom, search nearby for an aromatic atom.
        for sc_idx in side_chain_indices:
            if side_chain_has_aromatic(sc_idx, visited={alpha_idx, sc_idx}, max_depth=6):
                aromatic_sidechain_found = True
                break
        if aromatic_sidechain_found:
            break

    if not aromatic_sidechain_found:
        return False, "Alpha–amino acid backbone found but no aromatic ring connected in the side chain"

    return True, "Contains integrated alpha–amino acid backbone with an aromatic side chain"

# Example usage:
if __name__ == "__main__":
    # Test one of the examples: D-histidine
    test_smiles = "N[C@H](Cc1c[nH]cn1)C(O)=O"
    result, reason = is_aromatic_amino_acid(test_smiles)
    print(result, reason)