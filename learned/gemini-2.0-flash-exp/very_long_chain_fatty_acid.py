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

    # Check for carboxylic acid group or its deprotonated form
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,-1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found"


    max_chain_length = 0
    for acid_match in acid_matches:
      carbonyl_carbon_idx = acid_match[0]
      
      def get_chain_length(current_atom_idx, visited_atoms):
          """Recursively calculates the longest chain length."""
          
          max_len = 0
          
          for neighbor in mol.GetAtomWithIdx(current_atom_idx).GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor.GetSymbol() == "C" and neighbor_idx not in visited_atoms:
                
                new_visited = set(visited_atoms)
                new_visited.add(neighbor_idx)
                max_len = max(max_len, 1 + get_chain_length(neighbor_idx, new_visited))

          return max_len

      # Find neighboring carbons and start recursion
      carbonyl_carbon_neighbors = mol.GetAtomWithIdx(carbonyl_carbon_idx).GetNeighbors()
      current_chain_length = 0
      for neighbor in carbonyl_carbon_neighbors:
        if neighbor.GetSymbol() == "C":
             current_chain_length = max(current_chain_length, get_chain_length(neighbor.GetIdx(), {neighbor.GetIdx()}))
            
      max_chain_length = max(max_chain_length, current_chain_length + 1) # +1 to account for the carbon in the COOH

    if max_chain_length > 27:
        return True, f"Ultra-long-chain fatty acid (C={max_chain_length})"
    elif max_chain_length > 22:
        return True, f"Very long-chain fatty acid (C={max_chain_length})"
    else:
        return False, f"Not a very long-chain fatty acid (C={max_chain_length})"