"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: hydroxy fatty acid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid carrying one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[OH]')
    carboxy_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not carboxy_matches:
        return False, "No carboxylic acid group found"

    # Get the carbon atom(s) of the carboxylic acid group(s)
    carboxy_carbons = [match[0] for match in carboxy_matches]

    # Analyze each carboxylic acid group to check for a long hydrocarbon chain
    for carboxy_carbon_idx in carboxy_carbons:
        # Use BFS to find the longest carbon chain connected to the carboxylic acid carbon
        visited = set()
        to_visit = [(carboxy_carbon_idx, 0)]  # Tuple of (atom_idx, chain_length)
        max_chain_length = 0

        while to_visit:
            current_idx, chain_length = to_visit.pop()
            if current_idx in visited:
                continue
            visited.add(current_idx)
            atom = mol.GetAtomWithIdx(current_idx)

            # Only consider carbon atoms
            if atom.GetAtomicNum() == 6 and current_idx != carboxy_carbon_idx:
                chain_length += 1
                if chain_length > max_chain_length:
                    max_chain_length = chain_length

            # Add neighboring atoms to visit
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:
                    to_visit.append((neighbor_idx, chain_length))

        # Check if the chain length is sufficient (e.g., at least 4 carbons)
        if max_chain_length < 4:
            continue  # Try next carboxylic acid group
        else:
            break  # Valid chain found
    else:
        return False, f"No sufficient carbon chain found (max length {max_chain_length} carbons)"

    # Exclude atoms in the carboxylic acid group for hydroxy search
    carboxy_atoms = set()
    for match in carboxy_matches:
        carboxy_atoms.update(match)

    # Check for hydroxy groups (-OH) not part of the carboxylic acid
    hydroxy_pattern = Chem.MolFromSmarts('[CX4;!$(C=O)][OX2H]')
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    hydroxy_indices = [match[1] for match in hydroxy_matches if match[1] not in carboxy_atoms]

    if not hydroxy_indices:
        return False, "No hydroxy substituents found on the hydrocarbon chain"

    # All checks passed
    return True, "Molecule is a hydroxy fatty acid"