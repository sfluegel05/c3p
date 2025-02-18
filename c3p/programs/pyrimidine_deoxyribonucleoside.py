"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
"""
Classifies: pyrimidine deoxyribonucleoside
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside consists of a pyrimidine base attached to a deoxyribose sugar via a β-N1-glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for pyrimidine base (six-membered ring with nitrogens at positions 1 and 3)
    pyrimidine_pattern = Chem.MolFromSmarts("n1cnccc1")
    if not mol.HasSubstructMatch(pyrimidine_pattern):
        return False, "No pyrimidine base found"

    # Check for deoxyribose sugar (five-membered ring with oxygen and hydroxyl groups at positions 3' and 5')
    deoxyribose_pattern = Chem.MolFromSmarts("C1[C@@H]([C@H](O)[C@@H](CO)O1)")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No deoxyribose sugar found"

    # Check for β-N1-glycosidic bond between sugar and base
    glycosidic_bond_pattern = Chem.MolFromSmarts("n1c[C@H]2O[C@@H](C[C@H]2O)CO")
    if not mol.HasSubstructMatch(glycosidic_bond_pattern):
        return False, "No β-N1-glycosidic bond between sugar and base"

    # Ensure the absence of hydroxyl group at 2' position (i.e., it's deoxyribose, not ribose)
    ribose_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O")
    if mol.HasSubstructMatch(ribose_pattern):
        return False, "Contains ribose sugar instead of deoxyribose"

    # Check that the molecule is a nucleoside (base + sugar, no phosphate groups)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    if mol.HasSubstructMatch(phosphate_pattern):
        return False, "Molecule is a nucleotide, contains phosphate group"

    return True, "Contains pyrimidine base attached to deoxyribose via β-N1-glycosidic bond"