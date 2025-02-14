"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is a fatty acid with one or more hydroxy substituents.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for at least one carboxylic acid group (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"
    
    # Check for at least one hydroxyl group (-OH)
    alcohol_pattern = Chem.MolFromSmarts("[OH]")
    if not mol.HasSubstructMatch(alcohol_pattern):
         return False, "No hydroxy group found"

    # Get the carbons of the carboxylic acid groups
    carboxyl_carbon_matches = mol.GetSubstructMatches(acid_pattern)
    carboxyl_carbons = [match[0] for match in carboxyl_carbon_matches]

    def find_longest_chain_dfs(start_atom, visited_atoms, current_chain):
        """Depth-first search to find the longest carbon chain."""
        
        max_len = len(current_chain)
        max_path = list(current_chain)

        neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(start_atom).GetNeighbors() if n.GetAtomicNum() == 6]

        for neighbor in neighbors:
          if neighbor not in visited_atoms:
            visited_atoms.add(neighbor)
            len_path, path = find_longest_chain_dfs(neighbor, visited_atoms, current_chain + [neighbor])
            if len_path > max_len:
              max_len = len_path
              max_path = list(path)
            visited_atoms.remove(neighbor)
        
        return max_len, max_path
        
    # Check for a long aliphatic carbon chain connected to the carboxylic acid
    long_chain_found = False
    all_chain_atoms = set()
    for carboxyl_carbon in carboxyl_carbons:
      chain_carbons = [n.GetIdx() for n in mol.GetAtomWithIdx(carboxyl_carbon).GetNeighbors() if n.GetAtomicNum() == 6]

      for carbon in chain_carbons:
        visited = {carboxyl_carbon, carbon}
        chain_length, chain = find_longest_chain_dfs(carbon, visited, [carboxyl_carbon, carbon])
        
        if chain_length > 4: #Chain length of 4 means 5 carbons, since we included the carbonyl carbon as the first.
            long_chain_found = True
            all_chain_atoms.update(chain)
    
    if not long_chain_found:
      return False, "No long carbon chain attached to carboxylic group found"


    # Check for at least ONE hydroxy group on the fatty acid chain
    hydroxy_on_chain = False
    for atom in mol.GetAtoms():
      if atom.GetAtomicNum() == 6:
        for neighbor in atom.GetNeighbors():
          if neighbor.GetSymbol() == 'O' and neighbor.GetTotalValence() == 2:
            # Check if the hydroxyl-bearing carbon is part of the long chain
            if atom.GetIdx() in all_chain_atoms:
                hydroxy_on_chain = True
                break
        if hydroxy_on_chain:
           break

    if not hydroxy_on_chain:
        return False, "Hydroxy group not on the fatty acid chain"
    
    # Check the number of carbons:
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 8:
        return False, "Too few carbons to be a fatty acid"


    return True, "Contains a carboxylic acid group, at least one hydroxy group on the fatty acid chain and a long carbon chain"