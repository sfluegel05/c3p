"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: glycosphingolipid
"""
from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide.
    The carbohydrate residue is attached by a glycosidic linkage to O-1 of the sphingoid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycosphingolipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the sphingoid base
    # Sphingoid base: long-chain amino alcohol with hydroxyl at C1 and amino at C2
    sphingoid_smarts = "[#8][CH][CH](N)[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3]"
    sphingoid_pattern = Chem.MolFromSmarts(sphingoid_smarts)

    # Define a general SMARTS pattern for the ceramide backbone
    # Ceramide: N-acylated sphingoid base
    ceramide_smarts = "[#8][CH][CH](NC(=O)[#6])[CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH2][CH3]"
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)

    # Check for sphingoid base or ceramide backbone
    has_sphingoid = mol.HasSubstructMatch(sphingoid_pattern)
    has_ceramide = mol.HasSubstructMatch(ceramide_pattern)

    if not (has_sphingoid or has_ceramide):
        return False, "No sphingoid base or ceramide backbone found"

    # Define SMARTS pattern for a sugar moiety (pyranose ring)
    sugar_ring_smarts = "O[C@H]1[C@H][C@H][C@H][C@H]O1"
    sugar_ring_pattern = Chem.MolFromSmarts(sugar_ring_smarts)

    # Check for sugar ring
    has_sugar_ring = mol.HasSubstructMatch(sugar_ring_pattern)

    if not has_sugar_ring:
        return False, "No sugar ring detected"

    # Define SMARTS pattern for glycosidic linkage to O-1 of sphingoid
    # Looking for oxygen connecting sphingoid C1 and sugar anomeric carbon
    glycosidic_smarts = "[#8]-[CH]-1~[O]-[#6]-[#6]-[#6]-[#6]-[#8]-1"
    glycosidic_pattern = Chem.MolFromSmarts(glycosidic_smarts)

    # Check for glycosidic linkage
    has_glycosidic_linkage = mol.HasSubstructMatch(glycosidic_pattern)

    if not has_glycosidic_linkage:
        return False, "No carbohydrate residue attached via glycosidic linkage to O-1"

    return True, "Contains sphingoid or ceramide backbone with sugar moiety attached via glycosidic linkage to O-1"