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

    # Define the steroid core with a branching point at C13
    steroid_core_pattern = Chem.MolFromSmarts("[C]1[C][C][C]2[C]1[C][C]([C])3[C]2[C][C][C]4[C]3[C][C]4")
    if not mol.HasSubstructMatch(steroid_core_pattern):
          return False, "No steroid core found."
    
    # Get the matches of the steroid core.
    matches = mol.GetSubstructMatches(steroid_core_pattern)

    # If no matches, not a steroid
    if not matches:
        return False, "No steroid core match found"

    # Iterate through all matches, if there is one match, it should be a steroid.
    for match in matches:
      # C13 is the 8th atom in the SMARTS
      c13_index_in_match = 7
      c13_atom_index = match[c13_index_in_match]
      c13_atom = mol.GetAtomWithIdx(c13_atom_index)
      
      # Get neighbors of C13
      neighbors_c13 = c13_atom.GetNeighbors()

      #find C17. This will be attached to C13 and C16. C16 is connected to C13
      c17_index = -1
      for neighbor1 in neighbors_c13:
        if neighbor1.GetAtomicNum() == 6:
          # Get neighbors of C16.
          c16_neighbors = neighbor1.GetNeighbors()
          for neighbor2 in c16_neighbors:
              # Find neighbor that is not C13
              if neighbor2.GetIdx() != c13_atom_index and neighbor2.GetAtomicNum() == 6:
                 c17_neighbors = neighbor2.GetNeighbors()
                 for neighbor3 in c17_neighbors:
                     if neighbor3.GetAtomicNum() == 8:
                       if neighbor2.GetChiralTag() == Chem.ChiralType.CHI_TETRAHEDRAL_CW: #check for alpha configuration
                           return True, "Contains a steroid core with a 17-alpha hydroxyl group"
      
    return False, "No steroid core with a 17-alpha-hydroxy group found."