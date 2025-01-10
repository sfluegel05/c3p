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
    
    # Define aromatic flavonoid unit pattern
    flavonoid_pattern = Chem.MolFromSmarts("c1c(O)cc(-c2ccc(O)cc2)c(=O)c(O)c1")
    
    # Look for at least two flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    if len(flavonoid_matches) < 2:
        return False, "Too few flavonoid units found"

    # Check for linkage between flavonoid units typically present in biflavonoids
    # Look for C-C or C-O-C linkage using wildcard and variable length to accommodate typical biflavonoid linkages
    linkage_pattern = Chem.MolFromSmarts("c1cc(-!@[C]-!@[a])ccc1")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)

    # Assume each flavonoid match should participate in at least one such complex linkage
    if len(linkage_matches) < 1:
        return False, f"No sufficient linkages between flavonoid units found; found {len(linkage_matches)}"
    
    # Check aromatic carbon count and presence of oxo groups
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 20:
        return False, "Too few carbon atoms; unlikely to be a biflavonoid"
    if o_count < 4:
        return False, "Too few oxygen atoms; unlikely to be a biflavonoid"
    
    return True, "Contains multiple flavonoid units linked in a manner consistent with biflavonoids"