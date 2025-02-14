"""
Classifies: CHEBI:18133 hexose
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hexose(smiles: str):
    """
    Determines if a molecule is a hexose based on its SMILES string.
    A hexose is a six-carbon monosaccharide which in its linear form contains either an aldehyde group at position 1 (aldohexose) or a ketone group at position 2 (ketohexose).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hexose, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for the hexose backbone.
    # This will be used to verify the main 6-carbon chain, before checking carbonyl groups.
    linear_backbone_pattern = Chem.MolFromSmarts("[CX4][CX4][CX4][CX4][CX4][CX4]") # 6-carbon chain

    if not mol.HasSubstructMatch(linear_backbone_pattern):
        return False, "Molecule does not contain a 6-carbon backbone."

    # Check for pyranose (6-membered ring with O) as part of the backbone
    pyranose_pattern = Chem.MolFromSmarts("O1[CX4][CX4][CX4][CX4][CX4]1") # 6-membered ring with one oxygen, and 5 carbons.
    pyranose_matches = mol.HasSubstructMatch(pyranose_pattern)

    # Check for furanose (5-membered ring with O) as part of the backbone
    furanose_pattern = Chem.MolFromSmarts("O1[CX4][CX4][CX4][CX4]1") # 5-membered ring with one oxygen, and 4 carbons
    furanose_matches = mol.HasSubstructMatch(furanose_pattern)

    # Check if aldehyde at C1 (only for linear hexoses, or open chain)
    if not pyranose_matches and not furanose_matches:
        aldehyde_pattern = Chem.MolFromSmarts("[CX4]([#6])[C](=[O])[H]")
        aldehyde_matches = mol.GetSubstructMatches(aldehyde_pattern)

         #Check if ketone at C2 (only for linear hexoses, or open chain)
        ketone_pattern = Chem.MolFromSmarts("[CX4]([#6])[C](=[O])[CX4]")
        ketone_matches = mol.GetSubstructMatches(ketone_pattern)

        if len(aldehyde_matches) == 0 and len(ketone_matches) == 0:
            return False, "Not a hexose. Molecule does not have an aldehyde at C1 or a ketone at C2"

    # Count carbons of the backbone, not all carbons
    c_count = len(mol.GetSubstructMatches(linear_backbone_pattern)[0])

    # Count carbons of the pyranose backbone, or furanose.
    if pyranose_matches:
        c_count = len(mol.GetSubstructMatches(pyranose_pattern)[0])
    if furanose_matches:
        c_count = len(mol.GetSubstructMatches(furanose_pattern)[0])

    if c_count != 6:
        return False, f"Molecule does not have 6 carbons in the backbone, has {c_count}"



    return True, "Molecule is a hexose (6 carbons with aldehyde at C1 or ketone at C2 or ring structure)"