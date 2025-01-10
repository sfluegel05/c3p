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
    
    # Flavonoid unit pattern (expanded for common substituted structures)
    flavonoid_unit_pattern = Chem.MolFromSmarts("c1c(O)cc2oc3ccc(O)cc3cc2o1") # adding specific hydroxyl patterns
    flavonoid_matches = mol.GetSubstructMatches(flavonoid_unit_pattern)

    if len(flavonoid_matches) < 2:
        return False, "Less than two flavonoid units found"
    
    # Improved linkage patterns including heteroatoms and extended aromatic linkages
    linkage_patterns = [
        Chem.MolFromSmarts("c-O-c"),  # connection via oxygen
        Chem.MolFromSmarts("c(=O)-c"),  # linkage through carboxyl group
        Chem.MolFromSmarts("c-c"),  # direct carbon linkage, aromatic
        Chem.MolFromSmarts("c1c(O)cc2c(c1)c(c(c(c2O)O)O)O") # extended aromatic linkage pattern
    ]
    
    # Check for biflavonoid linking with more diverse patterns
    for linkage_pattern in linkage_patterns:
        if mol.HasSubstructMatch(linkage_pattern):
            return True, "Contains flavonoid units linked by identified biflavonoid patterns"

    return False, "Failed to identify expected biflavonoid linkage between flavonoid units"