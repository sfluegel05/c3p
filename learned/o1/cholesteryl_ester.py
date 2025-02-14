"""
Classifies: CHEBI:17002 cholesteryl ester
"""
"""
Classifies: cholesteryl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cholesteryl_ester(smiles: str):
    """
    Determines if a molecule is a cholesteryl ester based on its SMILES string.
    A cholesteryl ester is a sterol ester obtained by formal condensation of the carboxy group
    of any carboxylic acid with the 3-hydroxy group of cholesterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cholesteryl ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the cholesteryl ester pattern
    # Cholesterol backbone with an ester linkage at the 3-hydroxy position
    cholesteryl_ester_smarts = """
    [#6;A]1-[#6;A]-[#6;A]-[#6;A]2(-[#1])-[#6;A](-[#6;A]3(-[#6;A](-[#6;A]4(-[#6;A]-[#6;A]-[#6;A]-[#6;A]-[#6;A]-4)-[#1])(-[#1])-[#6;A]-3-[#1])(-[#1])-[#6;A]-2-[#6;A]-1)-[#8]-[#6](=O)-[#6]
    """
    cholesteryl_ester_pattern = Chem.MolFromSmarts(cholesteryl_ester_smarts)
    if cholesteryl_ester_pattern is None:
        return False, "Invalid cholesteryl ester SMARTS pattern"

    # Check for substructure match
    if not mol.HasSubstructMatch(cholesteryl_ester_pattern):
        return False, "Cholesteryl ester pattern not found"

    return True, "Molecule is a cholesteryl ester (cholesterol esterified at 3-hydroxy position)"