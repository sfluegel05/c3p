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
    2. Hydroxyl group at position 7 (A ring)
    3. Substituent at position 3 (B ring attachment)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define isoflavone core with 7-OH using precise SMARTS:
    # [A ring]-OH at position 7 (ortho to pyran oxygen)
    # Connected to chromen-4-one system (carbonyl at position 4)
    # B ring substituent at position 3 (C3 in pyran)
    iso_core = Chem.MolFromSmarts("[OH]c1c([CH]2[CH]=[CH]c3ccccc23)ccc2c1c(=O)ccoc2")
    
    # Alternative pattern accounting for possible conjugation variations
    alt_iso = Chem.MolFromSmarts("[OH]c1cc(O)c2c(c1)oc(-[!H0])cc(=O)c2")
    
    # Validate core structure first
    if mol.HasSubstructMatch(iso_core) or mol.HasSubstructMatch(alt_iso):
        # Ensure no esterification of the 7-OH
        ester_check = Chem.MolFromSmarts("[OH]c1ccc2c(c1)oc(-[!H0])cc(=O)c2-[OX2]")
        if not mol.HasSubstructMatch(ester_check):
            return True, "7-hydroxy group on isoflavone core"
    
    return False, "Does not match 7-hydroxyisoflavone structural requirements"