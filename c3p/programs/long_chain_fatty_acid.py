"""
Classifies: CHEBI:15904 long-chain fatty acid
"""
"""
Classifies: CHEBI:15904 long-chain fatty acid
A long-chain fatty acid is defined as a fatty acid with a chain length ranging from C13 to C22.
"""
from rdkit import Chem

def is_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a long-chain fatty acid based on its SMILES string.
    A long-chain fatty acid has a carboxylic acid group and a carbon chain length between 13 and 22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty acid, False otherwise
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

    # Function to calculate the longest carbon chain
    def calculate_longest_carbon_chain(mol):
        longest_chain = 0
        for atom in mol.GetAtoms():
            if atom.GetAtomicNum() == 6:  # Carbon atom
                visited = set()
                stack = [(atom, 1)]
                while stack:
                    current_atom, chain_length = stack.pop()
                    visited.add(current_atom.GetIdx())
                    if chain_length > longest_chain:
                        longest_chain = chain_length
                    for neighbor in current_atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                            stack.append((neighbor, chain_length + 1))
        return longest_chain

    # Calculate the longest carbon chain length
    longest_chain = calculate_longest_carbon_chain(mol)
    if longest_chain < 13 or longest_chain > 22:
        return False, f"Chain length is {longest_chain}, must be between 13 and 22"

    # Check for excessive branching (optional, but fatty acids are typically linear)
    n_branches = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetDegree() > 2)
    if n_branches > 2:  # Allow minimal branching
        return False, "Excessive branching detected"

    return True, "Contains a carboxylic acid group and a carbon chain length between 13 and 22"