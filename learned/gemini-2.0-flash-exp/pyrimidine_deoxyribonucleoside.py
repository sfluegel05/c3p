"""
Classifies: CHEBI:68472 pyrimidine deoxyribonucleoside
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_pyrimidine_deoxyribonucleoside(smiles: str):
    """
    Determines if a molecule is a pyrimidine deoxyribonucleoside based on its SMILES string.
    A pyrimidine deoxyribonucleoside contains a pyrimidine base linked to a 2-deoxyribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a pyrimidine deoxyribonucleoside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Generic pyrimidine base pattern. This pattern is more generic than the previous one
    pyrimidine_base_pattern1 = Chem.MolFromSmarts("[n]1ccccc1")
    pyrimidine_base_pattern2 = Chem.MolFromSmarts("[n]1c[c]nc[c]1")
    if not (mol.HasSubstructMatch(pyrimidine_base_pattern1) or mol.HasSubstructMatch(pyrimidine_base_pattern2)):
       return False, "No pyrimidine base detected"

    # Deoxyribose sugar, more robust pattern
    deoxyribose_pattern = Chem.MolFromSmarts("[C]1[C][C](O)[C](CO)[C](O)1")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
         return False, "No deoxyribose sugar detected"

    # Check if there is ribose sugar instead of deoxyribose
    ribose_pattern = Chem.MolFromSmarts("[C]1[C](O)[C]([C](O)CO)[C]1")
    if mol.HasSubstructMatch(ribose_pattern):
        return False, "Ribose detected, not deoxyribose"

    # Glycosidic bond check: Implied by connection between pyrimidine and deoxyribose pattern
    #  no need for the glycosidic_bond_pattern as the previous tests ensure connection

    return True, "Pyrimidine deoxyribonucleoside detected"