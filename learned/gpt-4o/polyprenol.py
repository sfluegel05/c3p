"""
Classifies: CHEBI:26199 polyprenol
"""
from rdkit import Chem

def is_polyprenol(smiles: str):
    """
    Determines if a molecule is a polyprenol based on its SMILES string.
    A polyprenol consists of more than one isoprene unit (H-[CH2C(Me)=CHCH2]nOH),
    with a terminal hydroxyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a polyprenol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the more specific isoprene unit pattern with branching
    isoprene_pattern = Chem.MolFromSmarts("[#6]=[#6]C([#6])[#6]")  # Represents C=CC(C)C
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than one"

    # Check for terminal hydroxyl group - must be connected at the end of the chain
    terminal_hydroxyl_pattern = Chem.MolFromSmarts("[OH]-[C]")  # Hydroxyl connected to carbon
    if not mol.HasSubstructMatch(terminal_hydroxyl_pattern):
        return False, "No terminal hydroxyl group (OH) was found in the structure"

    # Ensure the terminal hydroxyl is appropriately positioned
    hydroxyl_match = mol.GetSubstructMatches(terminal_hydroxyl_pattern)
    for match in hydroxyl_match:
        hydroxyl_atom = mol.GetAtomWithIdx(match[0])
        connected_carbon = hydroxyl_atom.GetNeighbors()[0]
        # Verify the attached carbon is part of an isoprene unit
        if any(connected_carbon.IsInRing() for match in isoprene_matches):
            return False, "Hydroxyl group not connected to terminal isoprene unit"

    return True, "Contains more than one isoprene unit with a terminal hydroxyl group, typical of polyprenols"