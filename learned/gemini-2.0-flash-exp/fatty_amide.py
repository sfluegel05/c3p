"""
Classifies: CHEBI:29348 fatty amide
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_fatty_amide(smiles: str):
    """
    Determines if a molecule is a fatty amide based on its SMILES string.
    A fatty amide is a monocarboxylic acid amide derived from a fatty acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty amide, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for the presence of an amide group (-C(=O)N-)
    amide_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide group found"
    
    # Check if the chain attached to the carbonyl carbon is not a carboxyl
    carbonyl_pattern = Chem.MolFromSmarts("[CX3](=[OX1])N")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)

    for match in carbonyl_matches:
      carbonyl_atom_index = match[0]
      carbonyl_atom = mol.GetAtomWithIdx(carbonyl_atom_index)
    
      # find all the carbon neighbours
      carbon_neighbors = [neighbor.GetIdx() for neighbor in carbonyl_atom.GetNeighbors() if mol.GetAtomWithIdx(neighbor).GetAtomicNum() == 6]
    
      if not carbon_neighbors:
         return False, "No carbon attached to carbonyl"
      
      found_long_chain = False
      for neighbor_idx in carbon_neighbors:
          neighbor = mol.GetAtomWithIdx(neighbor_idx)
          
          #If the neighbor is part of a carbonyl, then it is not part of a fatty acid chain
          if neighbor.GetSymbol() == 'C' and neighbor.GetTotalValence() == 3 and any(n.GetSymbol() == 'O' for n in neighbor.GetNeighbors()):
              continue
            
          # Check for long carbon chain (at least 4 carbons)
          path = Chem.GetShortestPath(mol, carbonyl_atom_index, neighbor_idx)
          
          if len(path) >= 4:
              # Count carbons in the attached chain.
              chain_carbon_count = 0
              for idx in path:
                if mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 :
                   chain_carbon_count += 1
              if chain_carbon_count >= 4:
                  found_long_chain = True
                  break
      
      if not found_long_chain:
          return False, "No long carbon chain (>=4) attached to the amide carbonyl"
    
    return True, "Molecule is a fatty amide"