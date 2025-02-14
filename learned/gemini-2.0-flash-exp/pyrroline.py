"""
Classifies: CHEBI:23763 pyrroline
"""
"""
Classifies: CHEBI:35419 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule is a pyrroline based on its SMILES string.
    A pyrroline is a 5-membered ring with one nitrogen and one double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrroline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for a 5-membered ring with one nitrogen and 1 double bond (anywhere in the ring)
    #This pattern will only match a 5 membered ring with a single N and a double bond.
    pyrroline_pattern = Chem.MolFromSmarts("[N]1[C]=[C]-[C]-[C]1")
    
    #Check for presence of the pattern
    if not mol.HasSubstructMatch(pyrroline_pattern):
      return False, "No pyrroline core structure found"
    
    #Ensure that it's not a pyrrole
    pyrrole_pattern = Chem.MolFromSmarts("[nX1][cX3]=[cX2][cX2]=[cX2]")
    if mol.HasSubstructMatch(pyrrole_pattern):
        return False, "It's a pyrrole not a pyrroline"
    
    return True, "Pyrroline structure detected"