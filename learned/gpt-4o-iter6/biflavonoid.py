"""
Classifies: CHEBI:50128 biflavonoid
"""
from rdkit import Chem

def is_biflavonoid(smiles: str):
    """
    Determines if a molecule is a biflavonoid based on its SMILES string.
    A biflavonoid is defined as a flavonoid oligomer obtained by the oxidative coupling 
    of at least two units of aryl-substituted benzopyran rings or its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a biflavonoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Improved SMARTS pattern for flavonoid unit (generic benzopyran substitute)
    flavonoid_unit_pattern = Chem.MolFromSmarts("c1cc2oc3ccccc3cc2o1")
    
    # Efficient way to find multiple flavonoid units
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_unit_pattern)
    if len(flavonoid_matches) < 2:
        return False, "Less than two flavonoid units found"
    
    # Potential linkage patterns between flavonoid units
    linkage_patterns = [
        Chem.MolFromSmarts("c-c"),
        Chem.MolFromSmarts("c-O-c"),
        Chem.MolFromSmarts("c-c-c"),
        Chem.MolFromSmarts("c(~O)~c") # aromatic carbon linking via oxygen
    ]
    
    # Check if there's evidence of a valid biflavonoid linkage in the molecule
    for linkage_pattern in linkage_patterns:
        if mol.HasSubstructMatch(linkage_pattern):
            return True, "Contains flavonoid units linked by expected biflavonoid linkages"

    return False, "Failed to identify expected biflavonoid linkage between flavonoid units"