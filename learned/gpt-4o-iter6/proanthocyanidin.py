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

    # Pattern for hydroxyflavan units - focusing on phenolic groups and flavan structure
    hydroxyflavan_pattern = Chem.MolFromSmarts("c1c(O)c(O)cc(O)c1-c2c(O)cc(O)c(O)c2O")  # Flavan unit with hydroxy groups
    hydroxyflavan_matches = mol.GetSubstructMatches(hydroxyflavan_pattern)
    if len(hydroxyflavan_matches) < 2:
        return False, "Less than two hydroxyflavan units found"
    
    # 4→8 linkage (a common interflavan bond pattern)
    linkage_4_8_pattern = Chem.MolFromSmarts("[C@H]([c]1[c]([OH])[c]([OH])[c]([OH])[c]([c]1[OH]))-[c]2[c]([OH])[c]([OH])[c]([OH])[c]([c]2[OH])")
    linkage_4_8_matches = mol.GetSubstructMatches(linkage_4_8_pattern)
    
    # 4→6 linkage (another interflavan bond pattern)
    linkage_4_6_pattern = Chem.MolFromSmarts("[C@H]([c]1[c]([OH])[c]([OH])[c]([OH])[c]([c]1[OH]))-[c]2[c]([OH])[c]([OH])[c]2")
    linkage_4_6_matches = mol.GetSubstructMatches(linkage_4_6_pattern)
    
    # Check if at least one linkage type is present
    if len(linkage_4_8_matches) == 0 and len(linkage_4_6_matches) == 0:
        return False, "No appropriate interflavan linkages found"
    
    return True, "Contains multiple hydroxyflavan units with appropriate linkages"