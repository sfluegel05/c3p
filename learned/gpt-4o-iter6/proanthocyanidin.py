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

    # Proanthocyanidins usually have multiple aromatic rings and hydroxyl groups
    # Look for flavan skeleton with hydroxy groups, typically C6-C3-C6 unit
    flavan_pattern = Chem.MolFromSmarts("c1cc(O)ccc1[C@H]2Oc3cc(O)cc(O)c3[C@@H]2")
    flavan_matches = mol.GetSubstructMatches(flavan_pattern)
    if len(flavan_matches) < 2:
        return False, "Less than two hydroxyflavan units found"

    # Check for plausible interflavan linkages, such as 4->8 linkages
    linkage_pattern = Chem.MolFromSmarts("[c,C]1[c,C][c,C]([C@H]2O[c,C][c,C][c,C][c,C]2)[c,C][c,C]1")
    linkage_matches = mol.GetSubstructMatches(linkage_pattern)
    if len(linkage_matches) < 1:
        return False, "No interflavan linkages found"

    return True, "Contains multiple hydroxyflavan units with appropriate linkages"