"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Bile Acid
"""
from rdkit import Chem

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a hydroxy-cholanic acid with a steroid backbone,
    hydroxyl groups, and a carboxylic acid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bile acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Remove stereochemistry from the molecule
    mol_2D = Chem.Mol(mol.ToBinary())
    Chem.RemoveStereochemistry(mol_2D)

    # Define the cholanic acid skeleton without stereochemistry
    cholanic_acid_smiles = "CC(CCC(=O)O)C1CCC2C3CCC4C(C)CCC4C3CCC12C"
    skeleton_pattern = Chem.MolFromSmiles(cholanic_acid_smiles)

    if skeleton_pattern is None:
        return False, "Error in skeleton pattern"

    # Remove stereochemistry from the skeleton pattern
    skeleton_pattern_2D = Chem.Mol(skeleton_pattern.ToBinary())
    Chem.RemoveStereochemistry(skeleton_pattern_2D)

    # Match the cholanic acid skeleton without considering stereochemistry
    if not mol_2D.HasSubstructMatch(skeleton_pattern_2D):
        return False, "Molecule does not match cholanic acid skeleton"

    # Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if num_hydroxyls == 0:
        return False, "No hydroxyl groups found"

    # Check for carboxylic acid group (COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    return True, "Molecule matches bile acid structural features"