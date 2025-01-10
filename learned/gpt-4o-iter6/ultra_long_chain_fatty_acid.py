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
    
    # Identify the starting carbon for chain analysis (usually the carbon bonded to the carboxyl carbon)
    carboxyl_carbon_idx = matches[0][0]
    carbon_near_carboxyl = None
    for neighbor in mol.GetAtomWithIdx(carboxyl_carbon_idx).GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
            carbon_near_carboxyl = neighbor.GetIdx()
            break
    
    if carbon_near_carboxyl is None:
        return False, "No starting carbon found for chain calculation"

    # To determine the longest chain from a starting carbon, using breadth-first search (BFS)
    def bfs_longest_chain(start_idx):
        queue = [(start_idx, 0)]
        visited = set()
        max_length = 0

        while queue:
            current_atom, current_length = queue.pop(0)
            if current_atom in visited:
                continue
            visited.add(current_atom)

            # Max chain length update
            if current_length > max_length:
                max_length = current_length
            
            for neighbor in mol.GetAtomWithIdx(current_atom).GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    queue.append((neighbor.GetIdx(), current_length + 1))

        return max_length
    
    # Get the longest chain starting from the chosen initial carbon
    longest_chain = bfs_longest_chain(carbon_near_carboxyl)

    # Criteria for classification as an ultra-long-chain fatty acid
    if longest_chain > 27:
        return True, f"Contains carboxylic acid group and chain length is C{longest_chain}, which is greater than C27"
    else:
        return False, f"Chain length is C{longest_chain}, which is not greater than C27"