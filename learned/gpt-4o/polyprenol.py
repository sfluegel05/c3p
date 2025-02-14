"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is composed of more than one isoprene unit (H-[CH2C(Me)=CHCH2]nOH),
    with a terminal hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced isoprene unit detection: CH2=C(Me)CH2
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) <= 1:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than one"

    # Check for primary hydroxyl group connectivity
    oxygens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8]
    hydroxyl_connected_carbon = None
    for oxygen in oxygens:
        if oxygen.GetDegree() == 1:  # Primary alcohol
            connected_atom = oxygen.GetNeighbors()[0]
            if connected_atom.GetAtomicNum() == 6:  # Check if connected to carbon
                hydroxyl_connected_carbon = connected_atom
                break
    if hydroxyl_connected_carbon is None:
        return False, "No terminal hydroxyl group found"

    # Ensure it's part of a valid isoprene chain
    if not hydroxyl_connected_carbon.HasSubstructMatch(isoprene_pattern):
        return False, "Hydroxyl group not connected to an isoprene unit"

    return True, "Contains more than one isoprene unit with a terminal hydroxyl group"