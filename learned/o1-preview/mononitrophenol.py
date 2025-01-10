"""
Classifies: CHEBI:39362 mononitrophenol
"""
from rdkit import Chem

def is_mononitrophenol(smiles: str):
    """
    Determines if a molecule is a mononitrophenol based on its SMILES string.
    A mononitrophenol is a phenol carrying a single nitro substituent at an unspecified position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a mononitrophenol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define nitro group pattern
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')

    # Find all nitro groups
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    num_nitro_groups = len(nitro_matches)
    if num_nitro_groups != 1:
        return False, f"Must have exactly one nitro group, found {num_nitro_groups}"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"

    # Find phenolic oxygen atoms (both protonated and deprotonated)
    phenol_pattern = Chem.MolFromSmarts('[c][O;H1,-1]')
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    phenol_o_indices = [match[1] for match in phenol_matches]

    if not phenol_o_indices:
        return False, "No phenolic OH group found"

    # For each ring, check if it contains exactly one phenolic oxygen and one nitro group
    for ring in atom_rings:
        ring_set = set(ring)
        hydroxyl_count = 0
        nitro_count = 0

        # Check for phenolic oxygen in ring
        for o_idx in phenol_o_indices:
            o_atom = mol.GetAtomWithIdx(o_idx)
            for neighbor in o_atom.GetNeighbors():
                if neighbor.GetIdx() in ring_set:
                    hydroxyl_count +=1
                    break  # Only count once per ring

        # Check for nitro group attached to ring
        for nitro_match in nitro_matches:
            n_idx = nitro_match[0]
            n_atom = mol.GetAtomWithIdx(n_idx)
            for neighbor in n_atom.GetNeighbors():
                if neighbor.GetIdx() in ring_set:
                    nitro_count +=1
                    break  # Only count once per ring

        # Check if ring has exactly one phenolic oxygen and one nitro group
        if hydroxyl_count == 1 and nitro_count == 1:
            return True, "Contains phenol ring with one nitro group attached"

    return False, "No ring with exactly one phenolic OH and one nitro group found"