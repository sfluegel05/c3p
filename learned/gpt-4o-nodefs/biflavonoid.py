"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    Biflavonoids generally contain two flavonoid moieties connected by covalent bonds.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # General structure pattern for a flavonoid
    flavonoid_pattern = Chem.MolFromSmarts("c1cc[cH]c2c1C(=O)CO2")  # Reduced flavonoid pattern for simplicity

    # Check for at least two occurrences of the basic flavonoid pattern
    matches = mol.GetSubstructMatches(flavonoid_pattern)
    if len(matches) < 2:
        return False, "Less than two flavonoid units found"

    # Check for connections typical of biflavonoids, (simple C-C or C-O between flavonoid units)
    c_o_linkage_pattern = Chem.MolFromSmarts("[c]~[o][c]")
    c_c_linkage_pattern = Chem.MolFromSmarts("[c][c][c]")

    # Verify the presence of covalent linkage between flavonoid units
    has_co_linkage = mol.HasSubstructMatch(c_o_linkage_pattern)
    has_cc_linkage = mol.HasSubstructMatch(c_c_linkage_pattern)

    if has_co_linkage or has_cc_linkage:
        return True, "Contains dual flavonoid units connected by typical linkage patterns"
    
    return False, "Molecule does not consist of connected flavonoid units typical of biflavonoids"