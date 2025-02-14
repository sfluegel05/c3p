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

    # Find all nitro groups
    nitro_pattern = Chem.MolFromSmarts('[N+](=O)[O-]')
    nitro_matches = mol.GetSubstructMatches(nitro_pattern)
    num_nitro_groups = len(nitro_matches)
    if num_nitro_groups != 1:
        return False, f"Must have exactly one nitro group, found {num_nitro_groups}"

    # Find all phenolic OH groups
    phenol_pattern = Chem.MolFromSmarts('[c][OH]')
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    if len(phenol_matches) == 0:
        return False, "No phenolic OH group found"

    # Get indices of nitro nitrogen atoms
    nitro_n_indices = [match[0] for match in nitro_matches]

    # Get indices of phenolic carbon atoms
    phenol_c_indices = [match[0] for match in phenol_matches]

    # Get ring information
    ri = mol.GetRingInfo()

    # For each phenolic OH group
    for phenol_match in phenol_matches:
        phenol_c = phenol_match[0]
        # Get rings that the phenolic carbon is in
        rings = [list(ring) for ring in ri.AtomRings() if phenol_c in ring]
        # Check if any ring contains both the phenolic carbon and the nitro nitrogen
        for ring in rings:
            ring_set = set(ring)
            # Check if nitro group is attached to the ring
            for n_idx in nitro_n_indices:
                # Get atom indices in nitro group (nitrogen and oxygens)
                nitro_atoms = [n_idx] + [nbr.GetIdx() for nbr in mol.GetAtomWithIdx(n_idx).GetNeighbors()]
                if ring_set.intersection(nitro_atoms):
                    return True, "Contains phenol ring with one nitro group attached"
    return False, "Phenolic ring with nitro group not found"