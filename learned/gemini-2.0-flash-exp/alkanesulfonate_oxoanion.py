"""
Classifies: CHEBI:134249 alkanesulfonate oxoanion
"""
from rdkit import Chem

def is_alkanesulfonate_oxoanion(smiles: str):
    """
    Determines if a molecule is an alkanesulfonate oxoanion based on its SMILES string.
    An alkanesulfonate oxoanion is characterized by a single sulfonate group (-S(=O)(=O)[O-])
    attached to an alkane chain (carbon chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkanesulfonate oxoanion, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfonate group (-S(=O)(=O)[O-])
    sulfonate_pattern = Chem.MolFromSmarts("[S](=[O])(=[O])[O-]")
    sulfonate_matches = mol.GetSubstructMatches(sulfonate_pattern)
    if not sulfonate_matches:
        return False, "No sulfonate group found"
    
    if len(sulfonate_matches) != 1:
        return False, f"Molecule contains {len(sulfonate_matches)} sulfonate groups, only 1 is allowed"
    
    sulfur_atom_idx = sulfonate_matches[0][0]
    sulfur_atom = mol.GetAtomWithIdx(sulfur_atom_idx)

    # Check for carbon directly attached to sulfonate sulfur
    carbon_found = False
    for neighbor in sulfur_atom.GetNeighbors():
       if neighbor.GetAtomicNum() == 6:
           carbon_found = True
           carbon_atom = neighbor
           break #exit the loop when a carbon is found
    if not carbon_found:
       return False, "No carbon directly attached to sulfonate sulfur"

    #Check for other sulfonate groups directly attached to that carbon
    for neighbor in carbon_atom.GetNeighbors():
      if neighbor.GetAtomicNum() == 16 and neighbor.GetIdx() != sulfur_atom_idx:
          return False, "Carbon attached to sulfonate contains another sulfonate"

    return True, "Molecule contains a single sulfonate group attached to a carbon chain"