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

    # Define a more comprehensive aromatic flavonoid substructure
    flavonoid_pattern = Chem.MolFromSmarts("c1cc(c(O)cc1)c2c(O)cc(O)c(C(=O)c2O)")
    
    # Look for at least two flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_pattern)
    if len(flavonoid_matches) < 2:
        return False, "Too few flavonoid units found"

    # Define linkage pattern for a variety of C-C and C-O-C linkages
    c_c_o_linkage_pattern = Chem.MolFromSmarts("c:c(:c):c-c:c(:o):c")
    c_o_c_linkage_pattern = Chem.MolFromSmarts("c:c(:o):c-O-c:c(:c):c")
    
    # Check for linkages
    c_c_linkages = len(mol.GetSubstructMatches(c_c_o_linkage_pattern))
    c_o_linkages = len(mol.GetSubstructMatches(c_o_c_linkage_pattern))
    if c_c_linkages + c_o_linkages == 0:
        return False, "No sufficient linkages between flavonoid units found"

    # Check for sufficient oxygen functionality
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, "Too few oxygens; unlikely to be a biflavonoid"
    
    return True, "Contains multiple flavonoid units linked consistent with biflavonoids"