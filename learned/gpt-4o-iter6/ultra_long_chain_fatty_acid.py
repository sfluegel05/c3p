"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid based on its SMILES string.
    Ultra-long-chain fatty acids have a chain length greater than C27 and a carboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (None, "Invalid SMILES string")

    # Ensure the molecule has a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    if not matches:
        return False, "No carboxylic acid group found"
    
    # Assuming the first match is the relevant carboxylic acid for analysis.
    carboxyl_idx = matches[0][0]  # First carbon of the carboxyl group
    
    # Function to recursively determine the longest carbon chain from the carboxyl carbon
    def get_longest_chain_length(start_atom_idx, visited):
        stack = [(start_atom_idx, 0)]
        max_length = 0
        while stack:
            current_atom_idx, current_length = stack.pop()
            if current_atom_idx in visited:
                continue
            visited.add(current_atom_idx)

            if current_length > max_length:
                max_length = current_length

            for neighbor in mol.GetAtomWithIdx(current_atom_idx).GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:  # Carbon atom
                    stack.append((neighbor.GetIdx(), current_length + 1))
        return max_length
    
    # Start searching from the carboxyl carbon's connected carbon (if available)
    start_idx = None
    for neighbor in mol.GetAtomWithIdx(carboxyl_idx).GetNeighbors():
        if neighbor.GetAtomicNum() == 6:  # Connected carbon atom
            start_idx = neighbor.GetIdx()
            break
    
    if start_idx is None:
        return False, "No starting carbon found for chain calculation"

    # Get the longest chain of contiguous carbon atoms
    longest_chain = get_longest_chain_length(start_idx, set())

    # Check if the chain length is greater than C27
    if longest_chain > 27:
        return True, f"Contains carboxylic acid group and chain length is C{longest_chain}, which is greater than C27"
    else:
        return False, f"Chain length is C{longest_chain}, which is not greater than C27"