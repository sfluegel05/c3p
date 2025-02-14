"""
Classifies: CHEBI:33856 aromatic amino acid
"""
"""
Classifies: aromatic amino acid
"""
from rdkit import Chem

def is_aromatic_amino_acid(smiles: str):
    """
    Determines if a molecule is an aromatic amino acid based on its SMILES string.
    An aromatic amino acid is an amino acid whose structure includes an aromatic ring in its side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic amino acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify alpha carbons: carbons connected to both an amino group and a carboxyl group
    alpha_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue  # Skip non-carbon atoms
        neighbor_elements = [nbr.GetAtomicNum() for nbr in atom.GetNeighbors()]
        if 7 in neighbor_elements and 6 in neighbor_elements:
            # Carbon connected to nitrogen (possible amino group) and carbon (possible carboxyl carbon)
            has_amino_group = False
            has_carboxyl_group = False
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 7:
                    # Nitrogen atom possibly part of amino group
                    if nbr.GetDegree() <= 3:
                        has_amino_group = True
                elif nbr.GetAtomicNum() == 6:
                    # Neighbor carbon atom possibly part of carboxyl group
                    oxygen_count = sum(1 for nnbr in nbr.GetNeighbors() if nnbr.GetAtomicNum() == 8)
                    if oxygen_count >= 2:
                        has_carboxyl_group = True
            if has_amino_group and has_carboxyl_group:
                alpha_carbons.append(atom)

    if not alpha_carbons:
        return False, "No amino acid backbone found"

    # For each alpha carbon, check the side chain
    for alpha_c in alpha_carbons:
        alpha_c_idx = alpha_c.GetIdx()
        # Identify amino group and carboxyl group atoms to exclude
        exclude_idxs = {alpha_c_idx}
        for nbr in alpha_c.GetNeighbors():
            if nbr.GetAtomicNum() == 7:
                # Amino group nitrogen
                exclude_idxs.add(nbr.GetIdx())
                # Exclude hydrogens attached to nitrogen
                for hnbr in nbr.GetNeighbors():
                    if hnbr.GetAtomicNum() == 1:
                        exclude_idxs.add(hnbr.GetIdx())
            elif nbr.GetAtomicNum() == 6:
                # Carboxyl carbon
                exclude_idxs.add(nbr.GetIdx())
                for nnbr in nbr.GetNeighbors():
                    if nnbr.GetAtomicNum() in [8,1]:
                        # Oxygens and hydrogens of carboxyl group
                        exclude_idxs.add(nnbr.GetIdx())
        # Traverse the side chain starting from the alpha carbon
        side_chain_atoms = set()
        visited = set(exclude_idxs)
        def dfs(atom_idx):
            for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                if nbr_idx not in visited:
                    visited.add(nbr_idx)
                    side_chain_atoms.add(nbr_idx)
                    dfs(nbr_idx)
        dfs(alpha_c_idx)

        if not side_chain_atoms:
            continue  # No side chain

        # Check if side chain contains aromatic ring(s)
        ring_info = mol.GetRingInfo()
        atom_rings = ring_info.AtomRings()
        if not atom_rings:
            continue  # No rings
        for ring in atom_rings:
            ring_set = set(ring)
            if ring_set & side_chain_atoms:
                # Side chain contains a ring
                # Check if all atoms in ring are aromatic
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    return True, "Aromatic ring found in side chain connected to amino acid backbone"
        # If no aromatic ring in side chain for this alpha carbon, continue to next
    return False, "Aromatic ring not found in side chain connected to amino acid backbone"

__metadata__ = {
    'chemical_class': {
        'name': 'aromatic amino acid',
        'definition': 'An amino acid whose structure includes an aromatic ring.'
    }
}