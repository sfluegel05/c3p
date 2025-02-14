"""
Classifies: CHEBI:2571 aliphatic alcohol
"""
from rdkit import Chem

def is_aliphatic_alcohol(smiles: str):
    """
    Determines if a molecule is an aliphatic alcohol based on its SMILES string.
    An aliphatic alcohol is an alcohol where at least one hydroxyl group is attached to a non-aromatic carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aliphatic alcohol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for the presence of at least one alcohol group (-OH) bound to an aliphatic carbon.
    alcohol_pattern = Chem.MolFromSmarts("[CX4H0-1][OX2H]")
    
    if not mol.HasSubstructMatch(alcohol_pattern):
      return False, "No aliphatic alcohol group found"

    #Verify it is an aliphatic carbon and not some other attachment.
    for match in mol.GetSubstructMatches(alcohol_pattern):
      for carbon_idx, oxygen_idx in zip(match[::2], match[1::2]): #The pattern has a form of C,O,C,O
          o_atom = mol.GetAtomWithIdx(oxygen_idx)
          c_atom = mol.GetAtomWithIdx(carbon_idx)
          
          if c_atom.GetIsAromatic():
            continue # Skip if its an aromatic carbon.
          
          #check for correct connectivity
          neighbors = o_atom.GetNeighbors()
          if len(neighbors) != 1:
            return False, "Oxygen not bound to a single carbon"
          if neighbors[0].GetIdx() != c_atom.GetIdx(): # should already be true given the SMARTS match, but lets double check
                return False, "Oxygen not bonded to a single carbon"
          if c_atom.GetAtomicNum() != 6:
              return False, "Alcohol not bound to a carbon"

          return True, "Molecule contains at least one alcohol group attached to a non-aromatic carbon."
    
    return False, "No aliphatic alcohol group found"