"""
Classifies: CHEBI:33855 arenecarbaldehyde
"""
"""
Classifies: arenecarbaldehyde
"""

from rdkit import Chem

def is_arenecarbaldehyde(smiles: str):
    """
    Determines if a molecule is an arenecarbaldehyde based on its SMILES string.
    An arenecarbaldehyde is any aldehyde in which the carbonyl group is attached to an aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an arenecarbaldehyde, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all aldehyde carbons: [C;H1](=O)
    aldehyde_pattern = Chem.MolFromSmarts('[C;H1](=O)')
    if aldehyde_pattern is None:
        return False, "Failed to create aldehyde pattern"

    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    n_aldehydes = len(aldehyde_matches)

    if n_aldehydes == 0:
        return False, "No aldehyde groups found in the molecule"

    # Initialize count of aldehydes attached to aromatic moieties
    n_aromatic_aldehydes = 0

    # For each aldehyde carbon, check if it is attached to an atom in an aromatic ring
    for match in aldehyde_matches:
        aldehyde_c_idx = match[0]
        aldehyde_c_atom = mol.GetAtomWithIdx(aldehyde_c_idx)

        # Get neighbor atoms excluding the oxygen
        neighbor_atoms = [nbr for nbr in aldehyde_c_atom.GetNeighbors() if nbr.GetAtomicNum() != 8]
        if not neighbor_atoms:
            continue  # No neighbors other than oxygen

        neighbor_atom = neighbor_atoms[0]

        # Check if neighbor atom is part of an aromatic ring
        if neighbor_atom.GetIsAromatic() and neighbor_atom.IsInRing():
            n_aromatic_aldehydes += 1

    if n_aromatic_aldehydes > 0:
        return True, f"Contains {n_aromatic_aldehydes} aldehyde group(s) attached to an aromatic moiety out of {n_aldehydes} aldehyde group(s) in the molecule"
    else:
        return False, f"No aldehyde group attached to an aromatic moiety found among {n_aldehydes} aldehyde group(s) in the molecule"