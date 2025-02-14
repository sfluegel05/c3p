"""
Classifies: CHEBI:143004 ultra-long-chain fatty acid
"""
from rdkit import Chem

def is_ultra_long_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is an ultra-long-chain fatty acid
    (chain length > C27) based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ultra-long-chain fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of a carboxylic acid group
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "Molecule does not contain a carboxylic acid group."

    # Find all carbon atoms
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]

    def dfs(current_node, visited):
      """Recursive Depth-First Search to find carbon chains"""
      max_path_length = 0
      
      for neighbor in mol.GetAtomWithIdx(current_node).GetNeighbors():
        if neighbor.GetAtomicNum() == 6:
          if neighbor.GetIdx() not in visited:
              path_length = dfs(neighbor.GetIdx(), visited | {neighbor.GetIdx()}) + 1
              max_path_length = max(max_path_length, path_length)

      return max_path_length
  
    max_chain_length = 0
    for start_node in carbon_atoms:
        chain_length = dfs(start_node, {start_node}) + 1 # +1 for the current atom
        max_chain_length = max(max_chain_length, chain_length)


    if max_chain_length > 27:
        return True, f"Longest carbon chain has {max_chain_length} carbons, which is > 27."
    else:
        return False, f"Longest carbon chain has {max_chain_length} carbons, which is <= 27."