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

    # Define amino acid backbone SMARTS pattern accounting for protonation states
    backbone_smarts = '[N;D2,D3,D4][C;H1,H2][C](=O)[O;H1,H0,-1]'
    pattern = Chem.MolFromSmarts(backbone_smarts)
    matches = mol.GetSubstructMatches(pattern)

    if not matches:
        return False, "No amino acid backbone found"

    # For each match, check for aromatic ring in side chain
    for match in matches:
        N_idx = match[0]
        C_alpha_idx = match[1]
        C_carboxyl_idx = match[2]
        O_carboxyl_idx = match[3]

        # Define backbone atom indices
        backbone_indices = set([N_idx, C_alpha_idx, C_carboxyl_idx, O_carboxyl_idx])

        # Include hydrogens attached to backbone atoms
        for idx in list(backbone_indices):
            atom = mol.GetAtomWithIdx(idx)
            for neighbor in atom.GetNeighbors():
                if neighbor.GetAtomicNum() == 1:
                    backbone_indices.add(neighbor.GetIdx())

        # Collect side chain atoms
        side_chain_indices = set()
        visited = set()
        stack = [C_alpha_idx]

        while stack:
            atom_idx = stack.pop()
            if atom_idx in visited or atom_idx in backbone_indices:
                continue
            visited.add(atom_idx)
            side_chain_indices.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    stack.append(neighbor_idx)

        if not side_chain_indices:
            continue  # No side chain atoms

        # Check if side chain contains an aromatic ring
        ring_info = mol.GetRingInfo()
        ring_atom_indices = ring_info.AtomRings()
        if not ring_atom_indices:
            continue  # No rings in molecule

        for ring in ring_atom_indices:
            ring_set = set(ring)
            if ring_set & side_chain_indices:
                # Side chain contains a ring
                # Check if all atoms in ring are aromatic
                if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                    return True, "Contains aromatic ring in side chain connected to amino acid backbone"

    # If no aromatic ring found in side chain for any backbone match
    return False, "Aromatic ring not found in side chain connected to amino acid backbone"

__metadata__ = {
    'chemical_class': {
        'name': 'aromatic amino acid',
        'definition': 'An amino acid whose structure includes an aromatic ring.'
    }
}