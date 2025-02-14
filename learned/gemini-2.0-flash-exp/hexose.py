"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide with either an aldehyde at position 1 or a ketone at position 2 in its linear form.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Basic check: 6 carbons and 6 oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count != 6:
        return False, f"Molecule does not have 6 carbons, has {c_count}"
    if o_count < 5 or o_count > 8: # Hexoses have between 5 and 6 oxygens, but substituted forms may have more
         return False, f"Molecule does not have between 5 and 8 oxygens, has {o_count}"

    # Check for linear form with aldehyde (aldose) at C1
    aldehyde_pattern = Chem.MolFromSmarts("[C](=[O])[C]") # carbon bonded to double bonded oxygen followed by another carbon.
    aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)


    # Check for linear form with ketone (ketose) at C2
    ketone_pattern = Chem.MolFromSmarts("[C](=[O])[C]") # carbon with double bonded oxygen followed by another carbon on each side.
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)

    #Check for pyranose (6-membered ring with O)
    pyranose_pattern = Chem.MolFromSmarts("O1[CX4][CX4][CX4][CX4][CX4]1") # 6-membered ring with one oxygen, and 5 carbons.
    pyranose_matches = mol.HasSubstructMatch(pyranose_pattern)
    #Check for furanose (5-membered ring with O)
    furanose_pattern = Chem.MolFromSmarts("O1[CX4][CX4][CX4][CX4]1") # 5-membered ring with one oxygen, and 4 carbons
    furanose_matches = mol.HasSubstructMatch(furanose_pattern)


    # Count Hydroxyl groups
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H1]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    # If not pyranose or furanose, check for aldose
    if not pyranose_matches and not furanose_matches and len(aldehyde_matches) == 0:
        if len(ketone_matches) == 0:
          return False, "Not a hexose (no aldehyde/ketone, and not a pyranose/furanose)"


    # If pyranose, must have 5 hydroxyls, or furanose must have 4 hydroxyls + 1 -CH2OH group
    if pyranose_matches and len(hydroxyl_matches) < 5:
        return False, "Pyranose ring does not have at least 5 hydroxyl groups"
    if furanose_matches and len(hydroxyl_matches) < 4:
         return False, "Furanose ring does not have at least 4 hydroxyl groups"


    return True, "Molecule is a hexose (6 carbons and at least one carbonyl or ring structure)"