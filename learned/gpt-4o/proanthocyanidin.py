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

    # Define catechin/epicatechin unit (hydroxyflavan) pattern
    hydroxyflavan_pattern = Chem.MolFromSmarts("c1cc(O)c2c(c1)C[C@@H](C2O)c3cc(O)cc(O)c3")
    
    # Check for presence of hydroxyflavan units
    flavan_matches = mol.GetSubstructMatches(hydroxyflavan_pattern)
    if len(flavan_matches) < 2:
        return False, f"Only {len(flavan_matches)} hydroxyflavan unit(s) found, need at least 2"

    # Check for 4→8 linkage pattern
    linkage_4_8_pattern = Chem.MolFromSmarts("[C@H]1(Oc2cc(O)c3c(oc2c(c3)O)C1)c4ccc(O)c(O)c4")
    if mol.HasSubstructMatch(linkage_4_8_pattern):
        return True, "Contains multiple hydroxyflavan units with 4→8 linkages typical of proanthocyanidins"
    
    # Check for 4→6 linkage pattern
    linkage_4_6_pattern = Chem.MolFromSmarts("[C@H]1(Oc2cc(O)c3c(oc2c(c3)O)C1)c4cc(O)cc(O)c4")
    if mol.HasSubstructMatch(linkage_4_6_pattern):
        return True, "Contains multiple hydroxyflavan units with 4→6 linkages typical of proanthocyanidins"

    return False, "Hydroxyflavan units present, but missing typical (4→8) or (4→6) linkages"