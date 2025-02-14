"""
Classifies: CHEBI:35746 fatty aldehyde
"""
"""
Classifies: fatty aldehyde (CHEBI:35581)
"""

from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is an aldehyde formally arising from reduction of the carboxylic acid group
    of its corresponding fatty acid, having a carbonyl group at one end of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define functional group SMARTS patterns
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H]=O")
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    ketone_pattern = Chem.MolFromSmarts("[CX3](=O)[#6]")

    # Check for exactly one aldehyde group
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)
    if len(aldehyde_matches) != 1:
        return False, f"Contains {len(aldehyde_matches)} aldehyde groups, expected exactly one."

    # Check for carboxylic acid groups
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) > 0:
        return False, "Contains carboxylic acid group(s), not a fatty aldehyde."

    # Check for ketone groups
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    # Exclude the aldehyde match from ketone matches
    ketone_atoms = set([atom_idx for match in ketone_matches for atom_idx in match])
    aldehyde_atoms = set([atom_idx for match in aldehyde_matches for atom_idx in match])
    ketone_only_atoms = ketone_atoms - aldehyde_atoms
    if len(ketone_only_atoms) > 0:
        return False, "Contains ketone group(s), not a fatty aldehyde."

    # Get the aldehyde carbon atom
    aldehyde_carbon_idx = aldehyde_matches[0][0]
    aldehyde_carbon = mol.GetAtomWithIdx(aldehyde_carbon_idx)

    # Check that aldehyde carbon is terminal (connected to only one carbon)
    neighbors = aldehyde_carbon.GetNeighbors()
    carbon_neighbor = None
    for atom in neighbors:
        if atom.GetAtomicNum() == 6:
            carbon_neighbor = atom
        else:
            return False, "Aldehyde carbon has non-carbon neighbor, not terminal."
    if carbon_neighbor is None:
        return False, "Aldehyde carbon has no carbon neighbor, not part of a chain."

    # Traverse the carbon chain starting from the neighbor of aldehyde carbon
    visited = set()
    def traverse_chain(atom):
        visited.add(atom.GetIdx())
        count = 1 if atom.GetAtomicNum() == 6 else 0
        for neighbor in atom.GetNeighbors():
            n_idx = neighbor.GetIdx()
            if n_idx not in visited and neighbor.GetAtomicNum() == 6:
                count += traverse_chain(neighbor)
        return count

    chain_length = traverse_chain(carbon_neighbor)
    # Include aldehyde carbon in chain length
    chain_length += 1

    # Check if chain length is at least 4 carbons
    if chain_length < 4:
        return False, f"Carbon chain too short ({chain_length} carbons), not a fatty aldehyde."

    return True, "Molecule is a fatty aldehyde with a terminal aldehyde group and appropriate carbon chain."