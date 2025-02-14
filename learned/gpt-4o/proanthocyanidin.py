"""
Classifies: CHEBI:26267 proanthocyanidin
"""
from rdkit import Chem

def is_proanthocyanidin(smiles: str):
    """
    Determines if a molecule is a proanthocyanidin based on its SMILES string.
    A proanthocyanidin is a flavonoid oligomer with two or more units of hydroxyflavans linked together.

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

    # Define hydroxyflavan (similar to catechin/epicatechin) unit pattern
    hydroxyflavan_pattern = Chem.MolFromSmarts("c1ccc2c(c1)C[C@@H]([C@H](O2)c3cc(O)cc(O)c3)O")
    
    # Check for presence of hydroxyflavan units
    flavan_matches = mol.GetSubstructMatches(hydroxyflavan_pattern)
    if len(flavan_matches) < 2:
        return False, f"Only {len(flavan_matches)} hydroxyflavan unit(s) found, need at least 2"

    # Check for 4->8 linkage pattern
    linkage_4_8_pattern = Chem.MolFromSmarts("[C@H]1(O[C@@H]2CC[c3cc(O)c(O)c(O)c32])Cc4ccc(O)cc4")
    if not mol.HasSubstructMatch(linkage_4_8_pattern):
        # Check for 4->6 linkage pattern
        linkage_4_6_pattern = Chem.MolFromSmarts("[C@H]1(O[C@@H]2CC[c3cc(O)c(O)c(O)c32])Cc4cc(O)ccc4")
        if not mol.HasSubstructMatch(linkage_4_6_pattern):
            return False, "Missing flavan linkages (4->8 or 4->6)"

    return True, "Contains multiple hydroxyflavan units with linkages typical of proanthocyanidins"