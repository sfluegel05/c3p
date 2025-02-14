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

    # Define the isoprene unit pattern, representing H-[CH2C(Me)=CH-CH2]nOH
    isoprene_pattern = Chem.MolFromSmarts("C(=C)C-C")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    if len(isoprene_matches) < 2:
        return False, f"Found {len(isoprene_matches)} isoprene units, need more than one"

    # Check for terminal hydroxyl group by iterating through oxygen atoms
    hydroxyl_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 8:  # Identifying oxygen atoms
            neighbors = atom.GetNeighbors()
            # Ensure it's a terminal hydroxyl group connected to a carbon
            if len(neighbors) == 1 and neighbors[0].GetAtomicNum() == 6:
                hydroxyl_found = True
                break

    if not hydroxyl_found:
        return False, "No terminal hydroxyl group (OH) was found in the structure"

    # Verify the hydroxyl group is directly linked to an isoprene unit
    connected_to_isoprene = any(mol.HasSubstructMatch(Chem.MolFromSmarts(f"{Chem.MolToSmarts(neighbors[0])}C(=C)C-C")) 
                                for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 for neighbors in [atom.GetNeighbors()])

    if not connected_to_isoprene:
        return False, "Hydroxyl group not correctly connected to terminal isoprene unit"

    return True, "Contains more than one isoprene unit with a terminal hydroxyl group, typical of polyprenols"