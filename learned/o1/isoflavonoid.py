"""
Classifies: CHEBI:50753 isoflavonoid
"""
"""
Classifies: isoflavonoid
"""

from rdkit import Chem

def is_isoflavonoid(smiles: str):
    """
    Determines if a molecule is an isoflavonoid based on its SMILES string.
    An isoflavonoid is any 1-benzopyran with an aryl substituent at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an isoflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the isoflavonoid core
    # Using atom map numbers to identify specific positions
    isoflavonoid_smarts = """
    [#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1
    -[#6]2:[#6](:[#6]:[#6]:[#6]:[#6]:2)
    =[#8]
    -[#6]3:[#6]:[#6]:[#6]:[#6]:[#6]:3
    """

    # Remove whitespace and newlines
    isoflavonoid_smarts = ''.join(isoflavonoid_smarts.split())

    # Create molecule from SMARTS
    isoflavonoid_pattern = Chem.MolFromSmarts(isoflavonoid_smarts)
    if isoflavonoid_pattern is None:
        return False, "Invalid isoflavonoid SMARTS pattern"

    # Search for isoflavonoid core
    matches = mol.GetSubstructMatches(isoflavonoid_pattern)
    if not matches:
        return False, "Does not contain isoflavonoid core structure"

    # For each match, check for aryl substituent at position 3
    for match in matches:
        # Map the atom indices from the SMARTS to the molecule
        mapping = dict(zip(isoflavonoid_pattern.GetSubstructMatch(isoflavonoid_pattern), match))

        # Atom indices in the SMARTS pattern:
        # Positions 1, 2, and 3 correspond to atoms where substituents can be attached
        # We need to check if there is an aryl group attached at position 3

        # Position 3 in the pattern corresponds to atom idx 5 (after removing hydrogens)
        position_3_idx = mapping.get(5)
        if position_3_idx is None:
            continue

        # Get the atom at position 3
        position_3_atom = mol.GetAtomWithIdx(position_3_idx)

        # Check neighbors for aryl substituent
        neighbors = position_3_atom.GetNeighbors()
        for neighbor in neighbors:
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in match:
                # Check if neighbor is part of an aryl group (aromatic ring)
                # Exclude hydrogen atoms
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIsAromatic():
                    # Check if the neighbor is part of an aromatic ring of size 6
                    ring_info = mol.GetRingInfo()
                    if ring_info.IsAtomInRingOfSize(neighbor_idx, 6):
                        return True, "Contains isoflavonoid core with aryl substituent at position 3"

    return False, "Does not have aryl substituent at position 3 of isoflavonoid core"