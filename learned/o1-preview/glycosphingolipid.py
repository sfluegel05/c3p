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

    # Define a general SMARTS pattern for the amide bond (N-acylation)
    amide_smarts = "[NX3][CX3](=O)[CX4,CX3]"
    amide_pattern = Chem.MolFromSmarts(amide_smarts)

    # Find all amide groups
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide (N-acylation) found"

    # Define a general SMARTS pattern for the sphingoid base (long-chain amino alcohol)
    # Nitrogen connected to a chain with at least one hydroxyl group
    sphingoid_smarts = "[NX3][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][CX4][OX2H1,OX1H0]"
    sphingoid_pattern = Chem.MolFromSmarts(sphingoid_smarts)

    # Check if sphingoid base is present
    has_sphingoid = mol.HasSubstructMatch(sphingoid_pattern)
    if not has_sphingoid:
        return False, "No sphingoid base (long-chain amino alcohol) found"

    # Define SMARTS pattern for sugar moiety (pyranose ring)
    sugar_ring_smarts = "[OX2H][C@@H]1[O,C][C@@H][C@@H][C@@H][O,C]1"
    sugar_ring_pattern = Chem.MolFromSmarts(sugar_ring_smarts)

    # Check for sugar ring
    has_sugar_ring = mol.HasSubstructMatch(sugar_ring_pattern)
    if not has_sugar_ring:
        return False, "No sugar moiety detected"

    # Define SMARTS pattern for glycosidic linkage to O-1 of sphingoid
    glycosidic_linkage_smarts = "[OX2H]-[C]"
    glycosidic_linkage_pattern = Chem.MolFromSmarts(glycosidic_linkage_smarts)

    # Check for glycosidic linkage between sphingoid base and sugar
    # Find oxygen connected to sphingoid carbon and sugar anomeric carbon
    found_glycosidic_linkage = False
    for match in mol.GetSubstructMatches(glycosidic_linkage_pattern):
        o_atom = mol.GetAtomWithIdx(match[0])
        neighbors = o_atom.GetNeighbors()
        has_sphingoid_carbon = False
        has_sugar_carbon = False
        for neighbor in neighbors:
            if neighbor.IsInRing():  # Sugar ring carbon
                has_sugar_carbon = True
            elif neighbor.GetAtomicNum() == 6:  # Carbon connected to sphingoid base
                has_sphingoid_carbon = True
        if has_sphingoid_carbon and has_sugar_carbon:
            found_glycosidic_linkage = True
            break

    if not found_glycosidic_linkage:
        return False, "No glycosidic linkage between sphingoid base and sugar moiety found"

    return True, "Contains sphingoid base with sugar moiety attached via glycosidic linkage to O-1"