"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: glycosphingolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

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

    # Define SMARTS pattern for sphingoid base (long-chain amino alcohol)
    # Sphingoid base: long hydrocarbon chain with an amino group and two hydroxyl groups
    sphingoid_smarts = "[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6]~[#6;R0]-[#6;R0]-[O]-[C]-[N]"
    sphingoid_pattern = Chem.MolFromSmarts(sphingoid_smarts)

    # Define SMARTS pattern for ceramide (sphingoid base with amide-linked fatty acid)
    ceramide_smarts = "[C;R0](=O)-[N]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]-[C]"
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)

    # Check for sphingoid base or ceramide backbone
    has_sphingoid = mol.HasSubstructMatch(sphingoid_pattern)
    has_ceramide = mol.HasSubstructMatch(ceramide_pattern)

    if not (has_sphingoid or has_ceramide):
        return False, "No sphingoid base or ceramide backbone found"

    # Define SMARTS pattern for glycosidic linkage to O-1 of sphingoid
    # Looking for an oxygen atom connected to a sugar ring (pyranose or furanose)
    glycosidic_smarts = "[C;R0]-[O]-[C;R]-[C;R]-[O;R]"
    glycosidic_pattern = Chem.MolFromSmarts(glycosidic_smarts)

    has_glycosidic_linkage = mol.HasSubstructMatch(glycosidic_pattern)

    if not has_glycosidic_linkage:
        return False, "No carbohydrate residue attached via glycosidic linkage to O-1"

    # Optionally, we can check for the sugar ring more specifically
    # Define SMARTS pattern for sugar rings (pyranose and furanose)
    sugar_ring_smarts = "[$([C@H]1O[C@@H](CO)[C@H](O)[C@@H](O)[C@H]1O),$([C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O)]"
    sugar_ring_pattern = Chem.MolFromSmarts(sugar_ring_smarts)

    has_sugar_ring = mol.HasSubstructMatch(sugar_ring_pattern)

    if not has_sugar_ring:
        return False, "No sugar ring detected"

    return True, "Contains sphingoid or ceramide backbone with sugar moiety attached via glycosidic linkage to O-1"