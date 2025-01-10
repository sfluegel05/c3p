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

    # Define SMARTS pattern for sphingoid base (long-chain amino diol)
    # Sphingoid base: HO-C-C(N)-C(OH)-[CH2]n-CH3, where n >=9
    sphingoid_smarts = "[O;H1][C;H1]-[C;H1]([N])[C;H]([O;H1])[C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H2][C;H3]"
    sphingoid_pattern = Chem.MolFromSmarts(sphingoid_smarts)
    
    # Define SMARTS pattern for ceramide (sphingoid base with amide-linked fatty acid)
    ceramide_smarts = "[C;H1](=O)-[N]-[C;H1]-[C;H1]([O;H1])[C;H]([O;H1])[C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;!R][C;H3]"
    ceramide_pattern = Chem.MolFromSmarts(ceramide_smarts)

    # Check for sphingoid base or ceramide backbone
    has_sphingoid = mol.HasSubstructMatch(sphingoid_pattern)
    has_ceramide = mol.HasSubstructMatch(ceramide_pattern)

    if not (has_sphingoid or has_ceramide):
        return False, "No sphingoid base or ceramide backbone found"

    # Define SMARTS pattern for glycosidic linkage to O-1 of sphingoid
    # Looking for an oxygen atom connected to a sugar ring (pyranose or furanose) and to C1 of sphingoid base
    # Oxygen atom connected to anomeric carbon of sugar [C;R] and to sphingoid O-1 [C;H1]
    glycosidic_smarts = "[C;R]-[O]-[C;H1]"
    glycosidic_pattern = Chem.MolFromSmarts(glycosidic_smarts)

    has_glycosidic_linkage = mol.HasSubstructMatch(glycosidic_pattern)

    if not has_glycosidic_linkage:
        return False, "No carbohydrate residue attached via glycosidic linkage to O-1"

    # Optionally, check for sugar ring more specifically
    # Define SMARTS pattern for pyranose ring
    sugar_ring_smarts = "C1[C,O][C,O][C,O][C,O][C,O]O1"
    sugar_ring_pattern = Chem.MolFromSmarts(sugar_ring_smarts)

    has_sugar_ring = mol.HasSubstructMatch(sugar_ring_pattern)

    if not has_sugar_ring:
        return False, "No sugar ring detected"

    return True, "Contains sphingoid or ceramide backbone with sugar moiety attached via glycosidic linkage to O-1"