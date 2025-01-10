"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
"""
Classifies: CHEBI:XXXXX ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    An ultra-long-chain fatty acid is defined as having a carbon chain length greater than C27.

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

    # Check for carboxylic acid group (-C(=O)OH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Find the longest carbon chain
    longest_chain_length = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:  # Carbon atom
            # Traverse the chain starting from this atom
            chain_length = 1
            visited = set()
            stack = [(atom, chain_length)]
            while stack:
                current_atom, current_length = stack.pop()
                visited.add(current_atom.GetIdx())
                if current_length > longest_chain_length:
                    longest_chain_length = current_length
                for neighbor in current_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                        stack.append((neighbor, current_length + 1))

    # Check if the longest chain is greater than C27
    if longest_chain_length > 27:
        return True, f"Longest carbon chain length is {longest_chain_length} (greater than C27)"
    else:
        return False, f"Longest carbon chain length is {longest_chain_length} (not greater than C27)"