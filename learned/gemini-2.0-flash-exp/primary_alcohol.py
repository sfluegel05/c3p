"""
Classifies: CHEBI:15734 primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_primary_alcohol(smiles: str):
    """
    Determines if a molecule is a primary alcohol based on its SMILES string.
    A primary alcohol is a compound in which a hydroxy group, -OH, is attached to a
    saturated carbon atom which has either three hydrogen atoms attached to it or
    only one other carbon atom and two hydrogen atoms attached to it.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS pattern to match any alcohol bond
    alcohol_pattern = Chem.MolFromSmarts("[CH0-3][OH]")
    
    matches = mol.GetSubstructMatches(alcohol_pattern)
    
    if not matches:
      return False, "No alcohol group found"
    
    for match in matches:
        carbon_index = match[0]
        carbon = mol.GetAtomWithIdx(carbon_index)
        
        carbon_neighbors = 0
        
        for neighbor in carbon.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                carbon_neighbors += 1
            
        hydrogens_count = carbon.GetTotalNumHs()
        
        if carbon_neighbors == 1:
          if hydrogens_count != 2:
            return False, "Not a primary alcohol: Has one carbon neighbor but not 2 hydrogens"
        elif carbon_neighbors == 0:
          if hydrogens_count != 3:
              return False, "Not a primary alcohol: Has no carbon neighbors but not 3 hydrogens"
        elif carbon_neighbors > 1:
          return False, "Not a primary alcohol: Has more than one carbon neighbor"
        
    return True, "All alcohol groups are primary"