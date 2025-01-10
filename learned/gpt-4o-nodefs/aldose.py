"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose typically contains an aldehyde group, with hydroxyl groups on carbon atoms,
    forming either open-chain or cyclic hemiacetal forms (pyranoses/furanoses).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is an aldose, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Aldehyde group pattern - check for linear aldose form
    aldehyde_pattern = Chem.MolFromSmarts("[CH](=O)")
    if mol.HasSubstructMatch(aldehyde_pattern):
        # Check for open chain with hydroxyls suitable for aldose
        oh_pattern = Chem.MolFromSmarts("[CX4H2][OH]")
        oh_matches = mol.GetSubstructMatches(oh_pattern)
        if len(oh_matches) >= 2:
            return True, "Open-chain form with aldehyde and sufficient hydroxyl groups"

    # Check for cyclic hemiacetal (pyranose/furanose) forms
    hemiacetal_pattern = Chem.MolFromSmarts("O[C@H]1[CH2][C@H][O][C@H]1")
    if mol.HasSubstructMatch(hemiacetal_pattern):
        return True, "Cyclic form consistent with aldose (hemiacetal pyranose/furanose)"

    return False, "Structure does not fit typical aldose characteristics"