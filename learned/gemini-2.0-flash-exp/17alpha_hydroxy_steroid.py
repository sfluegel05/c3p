"""
Classifies: CHEBI:35342 17alpha-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_17alpha_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is a 17alpha-hydroxy steroid based on its SMILES string.
    A 17alpha-hydroxy steroid has a hydroxyl group at carbon 17 with alpha stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 17alpha-hydroxy steroid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a generic steroid core pattern based on fused rings
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C]2[C]3[C]1[C][C]([C])([C])2[C][C]3") # Modified core pattern
    if not mol.HasSubstructMatch(steroid_core_pattern):
        return False, "No steroid core found."

    # Get the matches of the steroid core.
    matches = mol.GetSubstructMatches(steroid_core_pattern)
    
    if not matches:
      return False, "No steroid core match found"

    for match in matches:
      # C13 is the 8th atom in the SMARTS
      c13_index_in_match = 7
      c13_atom_index = match[c13_index_in_match]
      c13_atom = mol.GetAtomWithIdx(c13_atom_index)
    
      # Get neighbors of C13
      neighbors_c13 = c13_atom.GetNeighbors()
      
      c17_atom = None
      for neighbor in neighbors_c13:
        if neighbor.GetAtomicNum() == 6:
           c17_atom = neighbor
           break
      
      if c17_atom is None:
        return False, "No C17 atom found connected to C13"
      
      c17_neighbors = c17_atom.GetNeighbors()

      for neighbor in c17_neighbors:
          if neighbor.GetAtomicNum() == 8: #Oxygen
              if neighbor.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW:  # Check for alpha configuration
                  return True, "Contains a steroid core with a 17-alpha hydroxyl group"
    
    return False, "No steroid core with a 17-alpha-hydroxy group found."