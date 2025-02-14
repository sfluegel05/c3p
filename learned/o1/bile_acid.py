"""
Classifies: CHEBI:3098 bile acid
"""
"""
Classifies: Bile Acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_bile_acid(smiles: str):
    """
    Determines if a molecule is a bile acid based on its SMILES string.
    A bile acid is a hydroxy-5β-cholanic acid with a steroid backbone,
    hydroxyl groups, 5β-configuration, and a carboxylic acid side chain.

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

    # Define steroid nucleus pattern (four fused rings with specific stereochemistry)
    steroid_nucleus_smarts = """
    [#6]1[#6][#6][#6][#6][#6]1
    [#6]2[#6][#6][#6][#6][#6]2
    [#6]3[#6][#6][#6][#6][#6]3
    [#6]4[#6][#6][#6][#6]4
    """

    steroid_nucleus_pattern = Chem.MolFromSmarts(''.join(steroid_nucleus_smarts.split()))
    if not mol.HasSubstructMatch(steroid_nucleus_pattern):
        return False, "Molecule does not have the steroid nucleus"

    # Check for 5β-configuration at ring junction (specific stereochemistry)
    # Define SMARTS for 5β-configuration
    configuration_pattern = Chem.MolFromSmarts("""
    [C@H]5([C@@H]([C@H]6[C@@H](CC5)CC[C@@H]7CCCC7)CC6)
    """)
    if not mol.HasSubstructMatch(configuration_pattern):
        return False, "Molecule does not have the 5β-configuration"

    # Check for carboxylic acid group (COOH) attached to side chain
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check for at least one hydroxyl group
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    num_hydroxyls = len(mol.GetSubstructMatches(hydroxyl_pattern))
    if num_hydroxyls == 0:
        return False, "No hydroxyl groups found"

    return True, "Molecule matches bile acid structural features"