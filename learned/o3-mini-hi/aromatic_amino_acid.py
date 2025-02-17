"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acid
Definition: An amino acid whose structure includes an aromatic ring.
A molecule qualifies if it contains an integrated (free) amino acid backbone,
i.e. if there is a carbon (the candidate “alpha” carbon) that is directly attached
to both a free (non‐amidated) amino group and a free carboxyl group, and if one of the
remaining substituents (the side chain) contains at least one aromatic atom.
This function uses a relaxed approach so as to capture both canonical aliphatic amino acids
(e.g. phenylalanine, tryptophan) and so–called aromatic amino acids such as anthranilic acid derivatives.
"""

from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines whether a molecule is an aromatic amino acid.
    The molecule is considered a free amino acid if there is at least one carbon atom (alpha carbon)
    that is connected to (a) a free (non-amidated) amine and (b) a free carboxyl group.
    Additionally, at least one non-backbone substituent of that carbon (the side chain)
    must contain an aromatic ring.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an aromatic amino acid, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens to help in counting
    mol = Chem.AddHs(mol)

    # Helper: check if a nitrogen attached to candidate alpha is a free amine.
    def is_free_amino(nitrogen, alpha_idx):
        # Must be nitrogen.
        if nitrogen.GetAtomicNum() != 7:
            return False
        # Count heavy‐atom neighbors (exclude hydrogens).
        heavy_neighbors = [nbr for nbr in nitrogen.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # One of these should be the alpha carbon.
        if alpha_idx not in [nbr.GetIdx() for nbr in heavy_neighbors]:
            return False
        # Allow for a free amine to be unsubstituted or mono‐substituted (e.g. N-methyl).
        # If there is a substituent other than the alpha carbon, check that it is not part of a carbonyl.
        extra = [nbr for nbr in heavy_neighbors if nbr.GetIdx() != alpha_idx]
        # Free amine: allow either no extra heavy atom or one extra that is not “amidating”
        # (i.e. not a carbonyl carbon). We check if an extra neighbor (if present) is bonded by a double bond
        # to an oxygen.
        if len(extra) > 1:
            return False
        if len(extra) == 1:
            extra_atom = extra[0]
            # If the extra heavy neighbor is carbon and is double-bonded to oxygen, it is likely an amide.
            for bond in nitrogen.GetBonds():
                if bond.GetOtherAtom(nitrogen).GetIdx() == extra_atom.GetIdx():
                    # Check bond type; if single, allow it
                    if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                        return False
        # Otherwise, if there are only hydrogens aside from the alpha connection, treat as free.
        return True

    # Helper: check if a carbon is a free carboxyl carbon (part of a free –COOH or –COO–),
    # when connected to candidate alpha (passed as alpha_idx).
    def is_free_carboxyl(carbon, alpha_idx):
        if carbon.GetAtomicNum() != 6:
            return False
        # Get heavy neighbors of this carbon.
        heavy_neighbors = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # Expect three heavy neighbors: one is the candidate alpha and two oxygens.
        if len(heavy_neighbors) != 3:
            return False
        if alpha_idx not in [nbr.GetIdx() for nbr in heavy_neighbors]:
            return False
        oxygens = [nbr for nbr in heavy_neighbors if nbr.GetAtomicNum() == 8]
        if len(oxygens) != 2:
            return False
        # One of the oxygens must be double-bonded to the carboxyl carbon (the carbonyl oxygen)
        # and the other is single-bonded (the hydroxyl or deprotonated oxygen).
        db_found = False
        sb_found = False
        for oxy in oxygens:
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), oxy.GetIdx())
            if bond is None:
                continue
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                db_found = True
            elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                sb_found = True
        if db_found and sb_found:
            return True
        return False

    # DFS to find if a substructure (starting at start_idx) contains an aromatic atom.
    def dfs_contains_aromatic(start_idx, backbone_set):
        visited = set()
        stack = [start_idx]
        while stack:
            current_idx = stack.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)
            if atom.GetIsAromatic():
                return True
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in visited and nbr.GetIdx() not in backbone_set:
                    stack.append(nbr.GetIdx())
        return False

    found_candidate_backbone = False

    # Iterate over carbon atoms as potential "alpha" carbons.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        alpha_idx = atom.GetIdx()
        # For an integrated amino acid backbone, the alpha carbon must be attached to a free amino group
        # and to a free carboxyl group.
        free_amino_found = None
        free_carboxyl_found = None
        # Iterate over heavy neighbors of candidate alpha.
        for nbr in atom.GetNeighbors():
            # Skip hydrogens.
            if nbr.GetAtomicNum() == 1:
                continue
            # Look for a free amine.
            if nbr.GetAtomicNum() == 7:
                if is_free_amino(nbr, alpha_idx):
                    free_amino_found = nbr
            # Look for a free carboxyl group.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != alpha_idx:
                if is_free_carboxyl(nbr, alpha_idx):
                    free_carboxyl_found = nbr
        # We require both free amino and carboxyl groups.
        if free_amino_found is None or free_carboxyl_found is None:
            continue

        # Mark that we have at least one alpha–amino acid backbone.
        found_candidate_backbone = True

        # Identify the side chain candidates: heavy neighbors of alpha that are not the free amine or carboxyl.
        backbone_set = {alpha_idx, free_amino_found.GetIdx(), free_carboxyl_found.GetIdx()}
        # Also add the oxygens from the carboxyl group.
        for oxy in free_carboxyl_found.GetNeighbors():
            if oxy.GetAtomicNum() == 8:
                backbone_set.add(oxy.GetIdx())
        side_chain_candidates = []
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                continue
            if nbr.GetIdx() not in backbone_set:
                side_chain_candidates.append(nbr.GetIdx())
        # If no distinct side chain exists (i.e. glycine), then skip.
        if not side_chain_candidates:
            continue

        # For each side chain candidate, check if any connected atom (excluding backbone) is aromatic.
        for sc in side_chain_candidates:
            if dfs_contains_aromatic(sc, backbone_set):
                return True, "Contains integrated amino acid backbone with an aromatic side chain"
        # If no side chain candidate gave an aromatic fragment, we report that backbone exists but side chain is not aromatic.
        # (We do not return immediately because we may find an alternative alpha carbon later.)
    if found_candidate_backbone:
        return False, "Alpha–amino acid backbone found but no aromatic ring in side chain"
    else:
        return False, "No integrated free alpha–amino acid backbone found"

# Example usage:
if __name__ == "__main__":
    # Test a couple of cases.
    test_smiles_list = [
        ("NC1=C(C=C(Cl)C=C1)C(O)=O", "2-amino-5-chlorobenzoic acid"),
        ("N[C@H](Cc1c[nH]cn1)C(O)=O", "D-histidine"),
        ("Nc1ccccc1C(O)=O", "anthranilic acid"),
        ("NC(Cc1ccc(O)cc1)C(O)=O", "tyrosine"),
        ("N[C@H](Cc1ccccc1)C(O)=O", "D-phenylalanine")
    ]
    for smiles, name in test_smiles_list:
        result, reason = is_aromatic_amino_acid(smiles)
        print(f"{name}: {result} ({reason})")