"""
Classifies: CHEBI:27283 very long-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_very_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a very long-chain fatty acid based on its SMILES string.
    A very long-chain fatty acid is a fatty acid with a carbon chain greater than C22.
    Ultra-long-chain fatty acids are a subset with chain length > C27.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a very long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Find the carbon atom of the carbonyl group
    carbonyl_carbon_match = mol.GetSubstructMatch(Chem.MolFromSmarts("C(=O)O"))
    if carbonyl_carbon_match:
      carbonyl_carbon_idx = carbonyl_carbon_match[0]
      
      # Get neighbors of the carbonyl carbon
      carbonyl_carbon_neighbors = mol.GetAtomWithIdx(carbonyl_carbon_idx).GetNeighbors()
      
      # Find the carbon atom adjacent to the carbonyl
      chain_start_carbon_idx = -1
      for neighbor in carbonyl_carbon_neighbors:
          if neighbor.GetSymbol() == "C":
              chain_start_carbon_idx = neighbor.GetIdx()
              break
      if chain_start_carbon_idx == -1:
            return False, "No carbon chain found"

      # Breadth-first search to find the longest chain of carbon atoms
      max_chain_length = 0
      queue = [(chain_start_carbon_idx, 1)]  # (atom index, current chain length)
      visited = set()

      while queue:
          current_idx, current_length = queue.pop(0)
          visited.add(current_idx)
          max_chain_length = max(max_chain_length, current_length)

          for neighbor in mol.GetAtomWithIdx(current_idx).GetNeighbors():
              neighbor_idx = neighbor.GetIdx()
              if neighbor.GetSymbol() == "C" and neighbor_idx not in visited:
                  queue.append((neighbor_idx, current_length + 1))

      # We subtract 1 from the total length because we didn't start at the carboxyl C.
      max_chain_length -= 1

      if max_chain_length > 27:
          return True, f"Ultra-long-chain fatty acid (C={max_chain_length})"
      elif max_chain_length > 22:
           return True, f"Very long-chain fatty acid (C={max_chain_length})"
      else:
         return False, f"Not a very long-chain fatty acid (C={max_chain_length})"
    else:
        return False, "No carboxyl group C atom found"