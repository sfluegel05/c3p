"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is typically a carboxylic acid with an aliphatic chain of 5 or fewer carbons.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a short-chain fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ensure molecule has a carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid functional group found"

    # Identify carbon atoms that are not part of aromatic rings or other complex structures
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    
    longest_chain_length = 0
    
    for carbon in carbons:
        # Perform a breadth-first search originating from the carbon atom connected to the carboxylic acid group
        queue = [(carbon.GetIdx(), 0)]
        visited = set()
        
        while queue:
            current_idx, length = queue.pop(0)
            
            if current_idx in visited:
                continue
            
            visited.add(current_idx)
            
            for neighbor in mol.GetAtomWithIdx(current_idx).GetNeighbors():
                if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() not in visited:
                    queue.append((neighbor.GetIdx(), length + 1))

            longest_chain_length = max(longest_chain_length, length)
    
    # A short-chain fatty acid has an aliphatic chain length of up to 5 excluding the carboxylic group carbon
    if longest_chain_length > 4:
        return False, f"Aliphatic chain too long ({longest_chain_length}), must be 4 or fewer excluding carboxyl carbon"

    return True, "Contains carboxylic acid group with a short aliphatic chain"