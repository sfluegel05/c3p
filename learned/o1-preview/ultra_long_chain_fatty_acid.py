"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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

    # Check for rings
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring structures"

    # Look for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid functional group found"

    # Check for exactly one carboxylic acid group
    if len(carboxy_matches) != 1:
        return False, f"Found {len(carboxy_matches)} carboxylic acid groups, expected exactly 1"

    # Get the carbon atom of the carboxylic acid group
    carboxylic_carbon_idx = carboxy_matches[0][0]
    carboxylic_carbon = mol.GetAtomWithIdx(carboxylic_carbon_idx)

    # Traverse the longest linear unbranched chain from the carboxylic carbon
    visited = set()
    max_chain_length = traverse_chain(carboxylic_carbon, visited)

    if max_chain_length > 27:
        return True, f"Contains linear hydrocarbon chain of {max_chain_length} carbons"
    else:
        return False, f"Linear hydrocarbon chain is {max_chain_length} carbons, which is not greater than 27"

def traverse_chain(atom, visited):
    """
    Recursively traverses the longest linear unbranched carbon chain starting from the given atom.

    Args:
        atom (rdkit.Chem.Atom): The starting atom.
        visited (set): Set of visited atom indices.

    Returns:
        int: Length of the longest linear chain from the starting atom.
    """
    if atom.GetAtomicNum() != 6:
        return 0  # Only consider carbon atoms

    idx = atom.GetIdx()
    if idx in visited:
        return 0
    visited.add(idx)

    # Exclude atoms that are not hydrogen, carbon, oxygen (for hydroxyl groups), or nitrogen (for unsaturation)
    allowed_atomic_nums = {1, 6, 7, 8}

    neighbor_lengths = []
    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx()
        if neighbor_idx in visited:
            continue
        if neighbor.GetAtomicNum() not in allowed_atomic_nums:
            return 0  # Disallow if other atoms are attached (e.g., branching)
        if neighbor.GetAtomicNum() == 6:
            bond = mol.GetBondBetweenAtoms(idx, neighbor_idx)
            if bond.GetBondType() not in [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE]:
                return 0  # Disallow triple bonds or unusual bond types
            length = traverse_chain(neighbor, visited.copy()) + 1  # Include current carbon
            neighbor_lengths.append(length)
        elif neighbor.GetAtomicNum() in [8, 7]:
            # Allow hydroxyl (-OH) or unsaturation (C=N) as substituents
            continue
        else:
            # Disallow other heteroatoms
            return 0

    if neighbor_lengths:
        return max(neighbor_lengths)
    else:
        return 1  # End of chain

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'ultra-long-chain fatty acid',
        'definition': 'Any very long-chain fatty acid which has a chain length greater than C27.',
        'parents': []
    }
}