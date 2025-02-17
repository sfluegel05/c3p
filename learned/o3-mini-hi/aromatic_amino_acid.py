"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acid
Definition: An amino acid whose structure includes an aromatic ring. In other words,
the molecule must contain an integrated (free) alpha–amino acid backbone
(with a free –NH2 attached to an alpha–carbon that also bears a free carboxyl group)
and that backbone must carry an aromatic side chain.
"""

from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    To qualify, the molecule must contain an alpha–amino acid backbone (an alpha–carbon
    that carries a free (non‐amidic) amino group and a free carboxyl group) and the side chain 
    (the substituent on the alpha–carbon that is not the amino or carboxyl group)
    must lead to an aromatic ring.
    
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
    
    # Add explicit hydrogens to count attached H's reliably.
    mol = Chem.AddHs(mol)
    
    # --- Step 1. Look for an alpha-amino acid backbone ---
    # Here we require:
    #  - an sp3 (non-aromatic) carbon (the alpha carbon) that has exactly one hydrogen.
    #  - that carbon is attached to a nitrogen (the free amine) and to a carboxyl group
    #    (a carbon bonded to a double-bonded oxygen and an –OH or –O^-).
    #
    # The SMARTS below is designed to recognize the alpha carbon.
    # Note: We require the alpha carbon to have exactly one hydrogen [H1].
    backbone_smarts = "[C;!a;H1](N)(C(=O)[O;H1,-])"
    aa_backbone = Chem.MolFromSmarts(backbone_smarts)
    matches = mol.GetSubstructMatches(aa_backbone)
    if not matches:
        return False, "Missing integrated alpha–amino acid backbone (free amino and acid groups on same carbon)"
    
    # --- Step 2. Check for aromatic side chain ---
    # For each alpha–carbon found in a matching backbone, get its neighbors.
    # Two of the neighbors are already assigned: the free amino (N) and the carboxyl carbon.
    # The remaining neighbor (if any) is the side chain.
    def side_chain_has_aromatic(start_idx, visited, max_depth=6):
        """
        Depth-first search from atom with index start_idx to see if any aromatic atom
        is encountered. Only traverses bonds; stops if max_depth reached.
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
    # In our SMARTS the order is assumed: (alpha-carbon, free amine, carboxyl carbon).
    for match in matches:
        alpha_idx, amine_idx, acid_idx = match
        alpha_atom = mol.GetAtomWithIdx(alpha_idx)
        neighbor_indices = [nbr.GetIdx() for nbr in alpha_atom.GetNeighbors()]
        # Exclude the amino and carboxyl atoms.
        side_chain_indices = [idx for idx in neighbor_indices if idx not in (amine_idx, acid_idx)]
        # We require that the residue has a side-chain (if not, it might be glycine which is not aromatic).
        if not side_chain_indices:
            continue
        # For each candidate side-chain atom, search in its vicinity (up to a few bonds) for an aromatic atom.
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
    # Test one of the true positive examples:
    test_smiles = "NC1=C(C=C(Cl)C=C1)C(O)=O"  # 2-amino-5-chlorobenzoic acid
    result, reason = is_aromatic_amino_acid(test_smiles)
    print(result, reason)