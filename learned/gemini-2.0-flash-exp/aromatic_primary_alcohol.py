"""
Classifies: CHEBI:33857 aromatic primary alcohol
"""
from rdkit import Chem

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

    # SMARTS pattern for an aromatic primary alcohol (aromatic-C-CH2-OH)
    aromatic_alcohol_pattern = Chem.MolFromSmarts("[cX3,cX2]~[C]~[OH]")

    # Find matches
    matches = mol.GetSubstructMatches(aromatic_alcohol_pattern)
    
    if not matches:
        return False, "No aromatic alcohol substructure found."

    # Verify that the C is primary alcohol (-CH2OH)
    for match in matches:
      for atom_index in match:
            atom = mol.GetAtomWithIdx(atom_index)
            if atom.GetAtomicNum() == 6: # Check carbon
                if atom.GetTotalDegree() == 3 and atom.GetTotalNumHs() == 2: # Primary C with 2 H
                    
                    for neighbor in atom.GetNeighbors():
                      if neighbor.GetAtomicNum() == 8:
                         # Check oxygen neighbor
                         if neighbor.GetTotalDegree() == 1:
                            return True, "Aromatic primary alcohol found."
    
    
    return False, "No aromatic primary alcohol found."