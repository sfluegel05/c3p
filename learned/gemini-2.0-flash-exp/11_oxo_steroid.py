"""
Classifies: CHEBI:47787 11-oxo steroid
"""
from rdkit import Chem

def is_11_oxo_steroid(smiles: str):
    """
    Determines if a molecule is an 11-oxo steroid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11-oxo steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the steroid core (fused 4 ring system, with no explicit bond specifications)
    steroid_core_smarts = "[C]1[C][C][C]2[C][C]([C]1)[C][C]3[C]([C]2)[C][C][C]4[C]3[C][C][C]4"
    steroid_core_pattern = Chem.MolFromSmarts(steroid_core_smarts)
    
    # Define SMARTS pattern for a carbonyl group (=O) attached to any carbon
    carbonyl_smarts = "[C]~[C](=O)"
    carbonyl_pattern = Chem.MolFromSmarts(carbonyl_smarts)

    if not mol.HasSubstructMatch(steroid_core_pattern):
         return False, "Not a steroid core."

    #Check if a carbonyl group is attached to any carbon in the steroid.
    matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not matches:
        return False, "No carbonyl group"
    
    #check for any carbonyl that is not in the steroid
    for match in matches:
        for atom_idx in match:
           atom = mol.GetAtomWithIdx(atom_idx)
           if atom.IsInRing():
              is_steroid_carbon = False
              for atom_idx_ring in mol.GetSubstructMatch(steroid_core_pattern):
                   if atom_idx == atom_idx_ring:
                      is_steroid_carbon= True
                      break
              if not is_steroid_carbon:
                continue
              
              # Get the neighbours of the carbonyl carbon
              carbon_neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
              # Get the atom of the carbonyl oxygen
              carbonyl_oxygen =  [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() == 'O' and atom_idx in [n.GetIdx() for n in a.GetNeighbors()]]

              
              if carbonyl_oxygen:
                return True, "Steroid with a carbonyl group at position 11"
              
    return False, "No 11-oxo steroid"