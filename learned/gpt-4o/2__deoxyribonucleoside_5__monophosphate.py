"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Pattern for 2'-deoxyribose sugar with hydroxyl at 3' position, no hydroxyl at 2'
    deoxyribose_pattern = Chem.MolFromSmarts("[C@@H]1([C@@H](O)[C@H]([C]([O-])O)O[C@@H]1CO)")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar backbone found"
    
    # Nucleobase patterns for adenine, guanine (purines) and thymine, cytosine, uracil (pyrimidines)
    purine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")
    pyrimidine_pattern = Chem.MolFromSmarts("n1c(=O)[nH]c(=O)c1")
    
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No common nucleobase found (purine or pyrimidine)"
    
    # Check for a phosphate group at 5'-location
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)([O-])O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No 5'-phosphate group found"

    return True, "Contains 2'-deoxyribose with 5'-phosphate and nucleobase"