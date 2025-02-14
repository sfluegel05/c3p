"""
Classifies: CHEBI:26666 short-chain fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_short_chain_fatty_acid(smiles: str):
    """
    Determines if a molecule is a short-chain fatty acid based on its SMILES string.
    A short-chain fatty acid is an aliphatic monocarboxylic acid with a chain length of less than C6
    (maximum of 5 carbons in the alkyl chain) and no non-hydrocarbon substituents other than
    the carboxylic acid group.

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

    # Check for non-hydrocarbon atoms (excluding H, C, O)
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num != 1 and atomic_num != 6 and atomic_num != 8:
            return False, "Non-hydrocarbon substituent found"

    # Check for carboxylic acid group with attached alpha carbon
    carboxylic_acid_pattern = Chem.MolFromSmarts("[C]-[C](=O)O")
    matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if not matches:
        return False, "No carboxylic acid group with attached alpha carbon found"


    # Get the carbon index connected to the carbonyl
    carbonyl_carbon_idx = matches[0][1] # index of the carbon in C(=O)O
    alpha_carbon_idx = matches[0][0] #index of carbon connected to the carbonyl

    
    # Check the carbon chain length
    chain_length = 0
    current_atom_idx = alpha_carbon_idx
    
    visited_atoms = set()
    
    while chain_length <= 5:
       
      visited_atoms.add(current_atom_idx)
      
      neighbors = []
      for neighbor in mol.GetAtomWithIdx(current_atom_idx).GetNeighbors():
          neighbor_idx = neighbor.GetIdx()
          if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited_atoms:
             neighbors.append(neighbor_idx)
      
      if len(neighbors) == 0:
          break
      elif len(neighbors) == 1:
           current_atom_idx = neighbors[0]
           chain_length += 1
      elif len(neighbors) > 1: # branch
           chain_length = 6 #stop the counting
           break; #this is not a simple alkyl chain.
    
    if chain_length > 5:
         return False, "More than 5 carbons in the alkyl chain"
   
    return True, "Short-chain fatty acid criteria met"