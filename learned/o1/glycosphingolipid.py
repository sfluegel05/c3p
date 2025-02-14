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
    A glycosphingolipid is a glycolipid that is a carbohydrate-containing derivative of a sphingoid or ceramide,
    where the carbohydrate is attached by a glycosidic linkage to O-1 of the sphingoid.

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

    # Define sphingoid base pattern: long-chain amino diol
    sphingoid_smarts = "[#6]-[CH](O)-[CH](NC(=O)?[*])[CH2]-[CH](O)[#6]"
    sphingoid_pattern = Chem.MolFromSmarts(sphingoid_smarts)
    if sphingoid_pattern is None:
        return False, "Error in defining sphingoid base pattern"

    sphingoid_matches = mol.GetSubstructMatches(sphingoid_pattern)
    if not sphingoid_matches:
        return False, "No sphingoid base found"

    # Define sugar (monosaccharide) pattern
    sugar_smarts = "[C@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1O"
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)
    if sugar_pattern is None:
        return False, "Error in defining sugar pattern"

    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if not sugar_matches:
        return False, "No sugar moiety found"

    # Check for glycosidic linkage between sugar and sphingoid base at O-1
    found_glycosidic_linkage = False
    for sph_match in sphingoid_matches:
        sph_O1_idx = sph_match[1]  # Index of O-1 in sphingoid base
        sph_O1_atom = mol.GetAtomWithIdx(sph_O1_idx)
        sph_O1_neighbors = [nbr.GetIdx() for nbr in sph_O1_atom.GetNeighbors()]

        for sugar_match in sugar_matches:
            sugar_anomeric_C_idx = sugar_match[0]  # Anomeric carbon in sugar
            sugar_C_atom = mol.GetAtomWithIdx(sugar_anomeric_C_idx)
            sugar_C_neighbors = [nbr.GetIdx() for nbr in sugar_C_atom.GetNeighbors()]

            # Check if O-1 of sphingoid is connected to anomeric carbon of sugar
            bond = mol.GetBondBetweenAtoms(sph_O1_idx, sugar_anomeric_C_idx)
            if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                found_glycosidic_linkage = True
                break
        if found_glycosidic_linkage:
            break

    if not found_glycosidic_linkage:
        return False, "No glycosidic linkage between sphingoid base and sugar at O-1"

    return True, "Molecule is a glycosphingolipid with sphingoid base and glycosidic linkage to sugar at O-1"