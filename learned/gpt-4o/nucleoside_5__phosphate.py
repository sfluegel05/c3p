"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    Identifies purine or pyrimidine base, ribosyl or deoxyribosyl sugar, and phosphate group at the 5' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Ribose and 2'-deoxyribose pattern
    sugar_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@H](CO)[C@@H](O)O1")
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](CO)[C@H](O)[C@@H](O)O1")

    has_sugar = mol.HasSubstructMatch(sugar_pattern) or mol.HasSubstructMatch(deoxyribose_pattern)
    
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"

    # Purine and Pyrimidine base pattern
    purine_pattern = Chem.MolFromSmarts("c1ncnc2ncnn2c1")
    pyrimidine_pattern = Chem.MolFromSmarts("c1cncnc1")

    has_base = mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)
    
    if not has_base:
        return False, "No purine or pyrimidine base found"

    # Correct phosphate group at the 5' position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    if not phosphate_matches:
        return False, "No phosphate group found at C-5'"

    return True, "Contains ribosyl or deoxyribosyl derivative of a purine/pyrimidine base with a 5'-phosphate group"