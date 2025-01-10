"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is a very long-chain fatty acid with a chain length greater than C27.
    Ultra-long-chain fatty acids are monocarboxylic acids with a linear hydrocarbon chain of more than
    27 carbons, possibly including unsaturation and hydroxyl groups, but without significant branching
    or cyclic structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid functional group found"

    # Check for exactly one carboxylic acid group
    if len(carboxy_matches) != 1:
        return False, f"Found {len(carboxy_matches)} carboxylic acid groups, expected exactly 1"

    # Get the carboxylic acid carbon atom index
    carboxy_c_idx = carboxy_matches[0][0]

    # Initialize variables
    max_chain_length = 0
    visited = set()

    def traverse_chain(atom_idx, prev_idx, current_length):
        atom = mol.GetAtomWithIdx(atom_idx)
        atom_symbol = atom.GetSymbol()

        # Allow carbon atoms and oxygen atoms bonded to carbon (e.g., hydroxyl groups)
        if atom_symbol == 'C' or (atom_symbol == 'O' and atom.GetDegree() == 1 and atom.GetNeighbors()[0].GetSymbol() == 'C'):
            pass
        else:
            return current_length

        visited.add(atom_idx)
        neighbors = atom.GetNeighbors()

        # Count only carbon atoms in the chain length
        if atom_symbol == 'C':
            current_length += 1

        # Exclude atoms with more than 2 heavy atom neighbors (branch points)
        heavy_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors) > 2:
            visited.remove(atom_idx)
            return current_length

        max_length = current_length
        for neighbor in neighbors:
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx == prev_idx:
                continue
            if neighbor_idx in visited:
                continue
            bond = mol.GetBondBetweenAtoms(atom_idx, neighbor_idx)
            # Only traverse single or double bonds
            if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE):
                continue
            # Recursively traverse neighbor atoms
            length = traverse_chain(neighbor_idx, atom_idx, current_length)
            if length > max_length:
                max_length = length

        visited.remove(atom_idx)
        return max_length

    # Start traversal from the carboxylic acid carbon atom
    chain_length = traverse_chain(carboxy_c_idx, -1, 0)

    if chain_length > 27:
        return True, f"Contains linear hydrocarbon chain of {chain_length} carbons"
    else:
        return False, f"Longest linear hydrocarbon chain is {chain_length} carbons, which is not greater than 27"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'ultra-long-chain fatty acid',
        'definition': 'Any very long-chain fatty acid which has a chain length greater than C27.',
        'parents': []
    }
}