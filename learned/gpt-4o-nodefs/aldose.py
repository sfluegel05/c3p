"""
Classifies: CHEBI:15693 aldose
"""
from rdkit import Chem

def is_aldose(smiles: str):
    """
    Determines if a molecule is an aldose based on its SMILES string.
    An aldose is typically characterized by an aldehyde group in the open chain
    and can form cyclic structures such as pyranoses or furanoses.

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

    # Check for open-chain form with aldehyde group
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[CH2][OH]")  # General aldehyde and hydroxyl match
    if mol.HasSubstructMatch(aldehyde_pattern):
        oh_pattern = Chem.MolFromSmarts("[CX4H1]([OH])[CX4][OH]")
        oh_matches = mol.GetSubstructMatches(oh_pattern)
        if len(oh_matches) >= 1:
            return True, "Open-chain form with aldehyde and hydroxyls"

    # Check for cyclic structures - pyranose or furanose forms
    # Check for pyranose form (six-membered ring with oxygens and hydroxyls)
    pyranose_pattern = Chem.MolFromSmarts("C1OC(CO)C(O)C1O")
    if mol.HasSubstructMatch(pyranose_pattern):
        return True, "Cyclic pyranose form detected consistent with aldose"

    # Check for furanose form (five-membered ring)
    furanose_pattern = Chem.MolFromSmarts("C1OC(CO)C1O")
    if mol.HasSubstructMatch(furanose_pattern):
        return True, "Cyclic furanose form detected consistent with aldose"

    return False, "Structure does not fit typical aldose characteristics"