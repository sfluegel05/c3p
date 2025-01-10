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

    # Define cholesteryl ester SMARTS pattern
    # Cholesterol backbone with ester linkage at 3-position
    cholesteryl_ester_smarts = """
    [C@H](OC(=O)[#6])[CH2][CH2][C@H]1[C@@H]([CH2][C@@H]2[C@H]([CH2][CH2]3[C@](C)([CH2][CH2][C@@H]3[C@@H]2CC=C1C)C)C)C
    """
    cholesteryl_ester_pattern = Chem.MolFromSmarts(cholesteryl_ester_smarts)
    if cholesteryl_ester_pattern is None:
        return False, "Invalid cholesteryl ester SMARTS pattern"

    # Check for substructure match with chirality
    if mol.HasSubstructMatch(cholesteryl_ester_pattern,useChirality=True):
        return True, "Molecule is a cholesteryl ester (cholesterol esterified at 3-hydroxy position)"
    else:
        return False, "Cholesteryl ester pattern not found"