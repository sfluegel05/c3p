"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_aromatic_primary_alcohol(smiles: str):
    """
    Determines if a molecule is an aromatic primary alcohol based on its SMILES string.
    An aromatic primary alcohol has a hydroxyl group attached to a -CH2- group which is itself
    attached to an aromatic ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aromatic primary alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for primary alcohol (-CH2-OH)
    primary_alcohol_pattern = Chem.MolFromSmarts("[CH2][OH]")
    
    # SMARTS pattern for aromatic carbon (lowercase c)
    aromatic_carbon_pattern = Chem.MolFromSmarts("c")
    
    # Find primary alcohol matches
    alcohol_matches = mol.GetSubstructMatches(primary_alcohol_pattern)

    if not alcohol_matches:
      return False, "No primary alcohol group found."

    # Iterate through each primary alcohol match and verify aromatic connection
    for match in alcohol_matches:
        alcohol_carbon_index = match[0] # index of the carbon in -CH2-OH
        
        # Check neighbors of the carbon of the primary alcohol (-CH2-OH)
        alcohol_carbon = mol.GetAtomWithIdx(alcohol_carbon_index)
        for neighbor in alcohol_carbon.GetNeighbors():
            if mol.HasSubstructMatch(Chem.MolFromSmarts(f"[{neighbor.GetSymbol()}]")) == False:
               continue
            if mol.HasSubstructMatch(aromatic_carbon_pattern, fromAtomIdx=neighbor.GetIdx()):
              return True, "Aromatic primary alcohol found."

    return False, "Primary alcohol not connected to an aromatic ring."