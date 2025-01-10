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

    # Generalized pattern for hydroxyflavan units (C6-C3-C6 with hydroxy groups)
    hydroxyflavan_pattern = Chem.MolFromSmarts("c1cc(O)c(cc1)[C@@H]2Oc3cc(O)c(O)cc3C2")
    hydroxyflavan_matches = mol.GetSubstructMatches(hydroxyflavan_pattern)
    if len(hydroxyflavan_matches) < 2:
        return False, "Less than two hydroxyflavan units found"

    # Patterns for common proanthocyanidin linkages
    # 4->8 linkage
    linkage_4_8_pattern = Chem.MolFromSmarts("[C@H]1(Oc2c(cccc2)[c,C]3[c,C][c,C][c,C][c,C]3)Oc4cc(O)ccc14")
    linkage_4_8_matches = mol.GetSubstructMatches(linkage_4_8_pattern)
    
    # 4->6 linkage 
    linkage_4_6_pattern = Chem.MolFromSmarts("[C@H]1(Oc2c(cc(O)c(O)c2)[c,C]3[c,C][c,C][c,C][c,C]3)Oc4cc(O)ccc14")
    linkage_4_6_matches = mol.GetSubstructMatches(linkage_4_6_pattern)
    
    if len(linkage_4_8_matches) + len(linkage_4_6_matches) < 1:
        return False, "No appropriate interflavan linkages found"

    return True, "Contains multiple hydroxyflavan units with appropriate linkages"