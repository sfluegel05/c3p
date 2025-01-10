"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer consisting of two or more hydroxyflavan units.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a proanthocyanidin, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Improved pattern for hydroxyflavan units - emphasizing aromatic structure and hydroxy groups
    hydroxyflavan_pattern = Chem.MolFromSmarts("c1cc(O)c(O)c(O)c1-c2c(O)cc(O)cc2O")
    hydroxyflavan_matches = mol.GetSubstructMatches(hydroxyflavan_pattern)
    if len(hydroxyflavan_matches) < 2:
        return False, "Less than two hydroxyflavan units found"

    # Improved patterns for common proanthocyanidin linkages
    # 4→8 linkage
    linkage_4_8_pattern = Chem.MolFromSmarts("[C@@H]1(O[C@H]2c3cc(O)cc(O)c3-c3cc(O)cc(O)c23)c2cc(O)cc(O)c12")
    linkage_4_8_matches = mol.GetSubstructMatches(linkage_4_8_pattern)
    
    # 4→6 linkage 
    linkage_4_6_pattern = Chem.MolFromSmarts("[C@@H]1(O[C@H]2c3cc(O)cc(O)c3-c3cc(O)c(O)cc23)c2cc(O)cc(O)c12")
    linkage_4_6_matches = mol.GetSubstructMatches(linkage_4_6_pattern)
    
    # Check if at least one linkage type exists
    if len(linkage_4_8_matches) + len(linkage_4_6_matches) < 1:
        return False, "No appropriate interflavan linkages found"

    return True, "Contains multiple hydroxyflavan units with appropriate linkages"