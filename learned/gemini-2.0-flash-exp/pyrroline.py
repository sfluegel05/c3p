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

    # Generalized pattern for a 5-membered ring with one N and one double bond
    # The nitrogen atom should be sp3 hybridized, hence the [NX3] pattern
    pyrroline_pattern = Chem.MolFromSmarts("[NX3]1~[#6]~[#6]~[#6]~[#6]1")
    
    if not mol.HasSubstructMatch(pyrroline_pattern):
        return False, "No pyrroline core structure found"
    
    #Ensure there is only one double bond within the ring
    double_bond_pattern = Chem.MolFromSmarts("[#6]=[#6]")
    matches = mol.GetSubstructMatches(double_bond_pattern)
    ring_matches = 0
    for match in matches:
      if mol.GetAtomWithIdx(match[0]).IsInRing() and mol.GetAtomWithIdx(match[1]).IsInRing():
        ring_matches += 1
    if ring_matches != 1:
      return False, "The pyrroline ring must have exactly 1 double bond."

    # Check that this is not a pyrrole, since a pyrrole contains 2 double bonds within the ring.
    pyrrole_pattern = Chem.MolFromSmarts("[nX2]1[cX3]=[cX2][cX2]=[cX2]1")
    if mol.HasSubstructMatch(pyrrole_pattern):
        return False, "It's a pyrrole not a pyrroline"

    return True, "Pyrroline structure detected"