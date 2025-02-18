"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: CHEBI:7-hydroxyisoflavones
"""
from rdkit import Chem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.
    Must have:
    1. Isoflavone core (benzopyran-4-one system)
    2. Hydroxyl group at position 7 (on the A ring adjacent to the pyrone oxygen)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define SMARTS pattern for 7-hydroxyisoflavone core:
    # - Benzopyran-4-one system (benzene fused to pyrone with ketone at position 4)
    # - Hydroxyl group on the benzene ring (A ring) adjacent to pyrone oxygen
    iso_pattern = Chem.MolFromSmarts("[OH]c1c2oc(=O)cc2ccc1")
    
    # Check if the pattern matches
    if mol.HasSubstructMatch(iso_pattern):
        return True, "7-hydroxy group on isoflavone core"
    
    # Alternative pattern allowing for possible substituents on the B ring
    alt_pattern = Chem.MolFromSmarts("[OH]c1ccc2c(c1)c(=O)oc-*")
    if mol.HasSubstructMatch(alt_pattern):
        return True, "7-hydroxy group on isoflavone core with B-ring substituent"
    
    return False, "Does not match 7-hydroxyisoflavone structural requirements"