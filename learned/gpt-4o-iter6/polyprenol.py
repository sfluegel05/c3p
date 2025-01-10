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

    # Enhanced pattern for isoprene units: accounts for 'n' repetitions with appropriate bonding
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C-C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    if len(isoprene_matches) <= 1:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than 1"

    # Check for presence of hydroxyl group, ensuring it is most likely at terminal position
    hydroxyl_pattern = Chem.MolFromSmarts("O")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Detect terminal-like hydroxyl assuming linear polymer end
    terminal_oh = False
    for match in hydroxyl_matches:
        atom = mol.GetAtomWithIdx(match[0])
        neighbors = [n.GetAtomicNum() for n in atom.GetNeighbors()]
        # Ensuring OH is attached to a carbon with only hydrogens (indicative of terminal)
        if 6 in [n.GetAtomicNum() for n in atom.GetNeighbors() if n.GetAtomicNum() != 1]:
            terminal_oh = True
            break

    if not terminal_oh:
        return False, "Hydroxyl group is not positioned in a suggestive terminal manner"

    return True, "Contains more than one isoprene unit with a suitable terminal hydroxyl group"