"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acid
Definition: An amino acid whose structure includes an aromatic ring.
A molecule qualifies if it contains an integrated (free) amino acid backbone,
i.e. if there is a candidate “alpha” carbon (an sp3 carbon) that is directly attached
to both a free (non‐amidated) amino group and a free carboxyl group, and 
if at least one of the side chain substituents (i.e. non‐backbone attachment of that carbon)
contains at least one aromatic atom.

To avoid misclassifying peptides (which have multiple amide bonds and tend to be larger),
we also require that the overall molecule have a relatively low number of heavy atoms.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines whether a molecule qualifies as a free aromatic amino acid.
    The algorithm is based on detecting one integrated free amino acid backbone.
    The candidate alpha carbon must be:
      - an sp3 carbon directly bonded to a free (non‐amidated) amine,
      - directly bonded to a free carboxyl group (C(=O)O, with one double bond and one single bond oxygen),
      - and it must have at least one remaining substituent (the “side chain”)
        in which a DFS search finds at least one aromatic atom.
    Additionally, the molecule must be reasonably small (in number of heavy atoms), 
    to exclude peptides.
    
    Args:
        smiles (str): SMILES string.
        
    Returns:
        bool: True if the molecule qualifies as a free aromatic amino acid, otherwise False.
        str: Explanation for the classification.
    """
    # Parse SMILES and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)

    # Use a heavy atom cutoff to help filter out peptides and larger molecules.
    heavy_atom_count = mol.GetNumAtoms() - sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 1)
    if heavy_atom_count > 30:
        return False, f"Molecule has {heavy_atom_count} heavy atoms; too large to be a free amino acid (likely a peptide or derivative)"

    # Define a helper to check for a free (non-amidated) amino group.
    def is_free_amino(nitrogen, alpha_idx):
        # Ensure the atom is nitrogen.
        if nitrogen.GetAtomicNum() != 7:
            return False
        # Get heavy neighbors (i.e. not hydrogens)
        heavy_nbrs = [nbr for nbr in nitrogen.GetNeighbors() if nbr.GetAtomicNum() > 1]
        # One heavy neighbor should be the candidate alpha carbon.
        if alpha_idx not in [nbr.GetIdx() for nbr in heavy_nbrs]:
            return False
        # Count the other heavy substituents.
        extra = [nbr for nbr in heavy_nbrs if nbr.GetIdx() != alpha_idx]
        # Allow no or one extra heavy neighbor.
        if len(extra) > 1:
            return False
        if len(extra) == 1:
            extra_atom = extra[0]
            # Check that the bond from nitrogen to the extra neighbor is a single bond.
            bond = nitrogen.GetBondBetweenAtoms(nitrogen.GetIdx(), extra_atom.GetIdx())
            if bond is None or bond.GetBondType() != rdchem.BondType.SINGLE:
                return False
        return True

    # Define a helper to check for a free carboxyl group.
    def is_free_carboxyl(carbon, alpha_idx):
        # Must be carbon. We expect the carboxyl carbon to be sp2.
        if carbon.GetAtomicNum() != 6 or carbon.GetHybridization() != rdchem.HybridizationType.SP2:
            return False
        # Get heavy neighbors (should be two oxygens and the alpha carbon).
        heavy_nbrs = [nbr for nbr in carbon.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_nbrs) != 3:
            return False
        # Check that one neighbor is the candidate alpha carbon.
        if alpha_idx not in [nbr.GetIdx() for nbr in heavy_nbrs]:
            return False
        # The remaining two should be oxygens.
        oxygens = [nbr for nbr in heavy_nbrs if nbr.GetAtomicNum() == 8 and nbr.GetIdx() != alpha_idx]
        if len(oxygens) != 2:
            return False
        # One oxygen must be bonded by a double bond (carbonyl) and the other by single bond (–OH or –O–).
        db_found = False
        sb_found = False
        for oxy in oxygens:
            bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), oxy.GetIdx())
            if not bond:
                continue
            if bond.GetBondType() == rdchem.BondType.DOUBLE:
                db_found = True
            elif bond.GetBondType() == rdchem.BondType.SINGLE:
                sb_found = True
        return db_found and sb_found

    # DFS search to check if a side chain fragment (excluding backbone atoms) contains any aromatic atom.
    def dfs_contains_aromatic(start_idx, backbone_set):
        visited = set()
        stack = [start_idx]
        while stack:
            idx = stack.pop()
            if idx in visited:
                continue
            visited.add(idx)
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetIsAromatic():
                return True
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() not in visited and nbr.GetIdx() not in backbone_set:
                    stack.append(nbr.GetIdx())
        return False

    candidate_found = False
    # Loop over candidate alpha carbons: sp3 carbons.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != rdchem.HybridizationType.SP3:
            continue
        alpha_idx = atom.GetIdx()
        
        free_amino_atom = None
        free_carboxyl_atom = None
        
        # Check neighbors of alpha for free amine and carboxyl.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                continue
            # Look for free amine.
            if nbr.GetAtomicNum() == 7 and is_free_amino(nbr, alpha_idx):
                free_amino_atom = nbr
            # Look for free carboxyl group.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != alpha_idx and is_free_carboxyl(nbr, alpha_idx):
                free_carboxyl_atom = nbr
        if free_amino_atom is None or free_carboxyl_atom is None:
            continue  # this alpha candidate does not have a complete free backbone
        
        # We found a candidate integrated amino acid backbone.
        candidate_found = True
        # Build the set of backbone atom indices: the alpha, the free amine and carboxyl (and carboxyl oxygens).
        backbone_set = {alpha_idx, free_amino_atom.GetIdx(), free_carboxyl_atom.GetIdx()}
        for nbr in free_carboxyl_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                backbone_set.add(nbr.GetIdx())
        # Identify side chain candidates as heavy neighbors of the alpha that are not in the backbone.
        side_chain_idxs = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() not in backbone_set]
        # If there is no side chain (e.g. glycine), then do not classify as aromatic.
        if not side_chain_idxs:
            continue
        
        # Check whether any side chain fragment (via DFS) contains an aromatic atom.
        for sc_idx in side_chain_idxs:
            if dfs_contains_aromatic(sc_idx, backbone_set):
                return True, "Contains integrated free amino acid backbone with an aromatic side chain"
        # If backbone found but side chain does not appear aromatic.
        return False, "Integrated free amino acid backbone found but side chain lacks an aromatic ring"
        
    if candidate_found:
        return False, "Integrated free amino acid backbone found but no aromatic ring in side chain"
    return False, "No integrated free alpha–amino acid backbone found"

# Example usage (you can run tests here):
if __name__ == "__main__":
    test_cases = [
        ("NC1=C(C=C(Cl)C=C1)C(O)=O", "2-amino-5-chlorobenzoic acid"),
        ("N[C@H](Cc1c[nH]cn1)C(O)=O", "D-histidine"),
        ("Nc1ccccc1C(O)=O", "anthranilic acid"),
        ("NC(Cc1ccc(O)cc1)C(O)=O", "tyrosine"),
        ("N[C@H](Cc1ccccc1)C(O)=O", "D-phenylalanine"),
        ("O=C(N[C@@H](CC1=CC=CC=C1)C(O)=O)[C@@H](NC(=O)[C@@H](N)C)CC2=CC=C(O)C=C2", "Ala-Tyr-Phe (dipeptide)")
    ]
    for smi, name in test_cases:
        res, reason = is_aromatic_amino_acid(smi)
        print(f"{name}: {res} -> {reason}")