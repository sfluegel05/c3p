"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is defined as having a carbon chain length greater than C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find the longest carbon chain starting from the carboxylic acid
    def find_longest_carbon_chain(mol):
        longest_chain = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon atom
                # Check if this carbon is part of a carboxylic acid group
                if any(neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 1 for neighbor in atom.GetNeighbors()):
                    visited = set()
                    stack = [(atom, 1)]  # (atom, current_chain_length)
                    while stack:
                        current_atom, chain_length = stack.pop()
                        visited.add(current_atom.GetIdx())
                        if chain_length > longest_chain:
                            longest_chain = chain_length
                        for neighbor in current_atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                                # Only follow the chain if it's not part of a ring or branch
                                if not neighbor.IsInRing() and neighbor.GetDegree() <= 2:
                                    stack.append((neighbor, chain_length + 1))
        return longest_chain

    longest_chain_length = find_longest_carbon_chain(mol)

    # Check if the longest chain length is greater than 22
    if longest_chain_length > 22:
        return True, f"Longest carbon chain length is {longest_chain_length}, which is greater than C22"
    else:
        return False, f"Longest carbon chain length is {longest_chain_length}, which is not greater than C22"