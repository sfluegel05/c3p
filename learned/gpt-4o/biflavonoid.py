"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is a flavonoid oligomer consisting of at least two flavonoid
    units linked together usually via C-C or C-O-C bonds.

    Args:
        smiles (str): The SMILES string of the molecule.

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more comprehensive aromatic flavonoid unit pattern
    flavonoid_pattern = Chem.MolFromSmarts("c1c(O)c(O)c2c(c1)c(=O)c(cc2)O")
    
    # Look for at least two flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    if len(flavonoid_matches) < 2:
        return False, "Too few flavonoid units found"

    # Improve detection for linkages: look for possible C-C or C-O-C patterns between flavonoid units
    linkage_pattern = Chem.MolFromSmarts("[cH]c(-c)c | [cH]c(-o)c")
    linkage_matches = [mol.GetSubstructMatches(Chem.MolFromSmarts(lp)) for lp in linkage_pattern]
    linkage_matches_count = sum(len(matches) for matches in linkage_matches) 
   
    # If we're looking for non-adjacent linkages, count linkage tests
    if linkage_matches_count < 1:
        return False, f"No sufficient linkages between flavonoid units found; found {linkage_matches_count}"
    
    # Check aromatic carbon count and presence of oxo or hydroxy groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbon atoms; unlikely to be a biflavonoid"
    if o_count < 4:
        return False, "Too few oxygen atoms; unlikely to be a biflavonoid"
    
    return True, "Contains multiple flavonoid units linked in a manner consistent with biflavonoids"