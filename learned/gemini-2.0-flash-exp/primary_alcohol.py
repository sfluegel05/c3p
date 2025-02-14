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
    
    # Modified SMARTS pattern that matches -C-OH,
    # where the carbon has to be saturated
    # and it has 0 or 1 carbon neighbor only
    primary_alcohol_pattern = Chem.MolFromSmarts("[CX4H0-3;!$(C=[C,N,O,S]);!$(C#[C,N])][OH]")
    carbon_pattern=Chem.MolFromSmarts("[CH0-3][OH]") #to identify the alcohol bond
    
    matches= mol.GetSubstructMatches(carbon_pattern)
    
    
    # Check if any matches are found
    if not matches:
        return False, "No alcohol group found"

    # Check each carbon to see if it's primary
    for match in matches:
      carbon_index = match[0]
      carbon = mol.GetAtomWithIdx(carbon_index)
      
      # count the number of carbon neighbors
      carbon_neighbors=0
      for neighbor in carbon.GetNeighbors():
        if neighbor.GetAtomicNum()==6:
          carbon_neighbors=carbon_neighbors+1

      if carbon_neighbors > 1:
            continue #not primary
      
      # now we verify that the matched carbon also has the pattern 
      if mol.HasSubstructMatch(primary_alcohol_pattern):
        return True, "Primary alcohol group found"
    return False, "No primary alcohol group found"