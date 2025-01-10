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

    # Identify amino acid backbone: [N]-[C@H]-[C](=O)-O
    amino_acid_pattern = Chem.MolFromSmarts("[N;D2]-[C;D3]([C;!H0,H1,H2,H3])[C](=O)-O")
    aa_matches = mol.GetSubstructMatches(amino_acid_pattern)
    if not aa_matches:
        return False, "No amino acid backbone found"

    # Assume first match is the amino acid backbone
    n_idx, alpha_c_idx, c_idx = aa_matches[0][0], aa_matches[0][1], aa_matches[0][2]

    # Find side chain attached to alpha carbon
    alpha_c = mol.GetAtomWithIdx(alpha_c_idx)
    side_chain_atom_indices = set()
    visited = set([n_idx, alpha_c_idx, c_idx])

    def dfs(atom_idx):
        for neighbor in mol.GetAtomWithIdx(atom_idx).GetNeighbors():
            nbr_idx = neighbor.GetIdx()
            if nbr_idx not in visited:
                visited.add(nbr_idx)
                side_chain_atom_indices.add(nbr_idx)
                dfs(nbr_idx)
    dfs(alpha_c_idx)

    # Remove backbone atoms from side chain atoms
    side_chain_atom_indices -= set([n_idx, alpha_c_idx, c_idx])

    if not side_chain_atom_indices:
        return False, "No side chain found"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"

    # Check if any ring in side chain is aromatic
    for ring in atom_rings:
        ring_set = set(ring)
        if ring_set & side_chain_atom_indices:
            # Check if all atoms in ring are aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                return True, "Aromatic ring found in side chain"
    return False, "Aromatic ring not found in side chain connected to amino acid backbone"

__metadata__ = {
    'chemical_class': {
        'name': 'aromatic amino acid',
        'definition': 'An amino acid whose structure includes an aromatic ring.'
    }
}