"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol is a prenol with more than one isoprene unit and a terminal hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for repeating isoprene units
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C-C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    if len(isoprene_matches) <= 1:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than 1"

    # Check for terminal hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(hydroxyl_pattern):
        return False, "Missing terminal hydroxyl group"

    # Check the position of the hydroxyl group to ensure it is at the end
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    terminal_oh = False
    for match in hydroxyl_matches:
        atom = mol.GetAtomWithIdx(match[0])
        neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
        # Ensure that OH is attached to a carbon (often the terminal position in polyprenols)
        if 6 in neighbors and len(neighbors) == 1:
            terminal_oh = True
            break

    if not terminal_oh:
        return False, "Hydroxyl group is not terminal"

    return True, "Contains more than one isoprene unit with a terminal hydroxyl group"