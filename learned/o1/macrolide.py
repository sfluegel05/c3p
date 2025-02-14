"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    A macrolide is a macrocyclic lactone with a ring of twelve or more members derived from a polyketide.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure ring information is computed
    mol.UpdatePropertyCache(strict=False)

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"

    # Define lactone pattern where both carbonyl carbon and ester oxygen are in a ring
    lactone_smarts = '[C;R](=O)[O;R]'
    lactone_pattern = Chem.MolFromSmarts(lactone_smarts)
    if lactone_pattern is None:
        return False, "Invalid lactone SMARTS pattern"

    # Find all matches of the lactone pattern
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)
    if not lactone_matches:
        return False, "No lactone functional groups found"

    # Iterate over rings of size >=12
    for ring_atoms in atom_rings:
        ring_size = len(ring_atoms)
        if ring_size >= 12:
            ring_atom_set = set(ring_atoms)
            # Check if any lactone match is fully within the ring
            for match in lactone_matches:
                if set(match).issubset(ring_atom_set):
                    return True, f"Contains macrocyclic lactone ring of size {ring_size} with lactone functionality"

    return False, "No macrocyclic lactone ring of size 12 or more with lactone functionality found"