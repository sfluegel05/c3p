"""
Classifies: CHEBI:28868 fatty acid anion
"""
"""
Classifies: fatty acid anion
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_anion(smiles: str):
    """
    Determines if a molecule is a fatty acid anion based on its SMILES string.
    A fatty acid anion is the conjugate base of a fatty acid, arising from deprotonation of the carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid anion, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Calculate total formal charge
    charge = sum([atom.GetFormalCharge() for atom in mol.GetAtoms()])
    if charge >= 0:
        return False, "Molecule does not have a net negative charge"

    # Look for carboxylate group (-C(=O)[O-])
    carboxylate_pattern = Chem.MolFromSmarts("[CX3](=O)[O-]")
    if not mol.HasSubstructMatch(carboxylate_pattern):
        return False, "No carboxylate group found"

    # Get the carbon atom of the carboxylate group
    carboxylate_matches = mol.GetSubstructMatches(carboxylate_pattern)
    carboxylate_carbons = [match[0] for match in carboxylate_matches]

    # Initialize maximum aliphatic chain length
    max_chain_length = 0

    # Check for aliphatic chain connected to carboxylate carbon
    for carboxylate_carbon in carboxylate_carbons:
        # Start from the carboxylate carbon and traverse the chain
        visited = set()
        stack = [(carboxylate_carbon, 0)]
        while stack:
            atom_idx, length = stack.pop()
            if atom_idx in visited:
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            # Skip carboxylate oxygens
            if atom.GetAtomicNum() == 8:
                continue
            # Update maximum chain length
            if length > max_chain_length:
                max_chain_length = length
            # Traverse neighbors
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                neighbor_atom = mol.GetAtomWithIdx(neighbor_idx)
                if neighbor_idx not in visited and neighbor_atom.GetAtomicNum() == 6 and not neighbor_atom.IsAromatic():
                    stack.append((neighbor_idx, length + 1))

    # Check that the longest aliphatic chain is at least 4 carbons long
    if max_chain_length < 4:
        return False, f"Aliphatic chain too short (length {max_chain_length}), need at least 4"

    # Check that the molecule is not aromatic
    if mol.HasSubstructMatch(Chem.MolFromSmarts("a")):
        return False, "Molecule contains aromatic rings"

    return True, "Molecule is a fatty acid anion with a carboxylate group and a long aliphatic chain"