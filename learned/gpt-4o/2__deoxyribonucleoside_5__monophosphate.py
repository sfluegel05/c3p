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

    # Pattern for 2'-deoxyribose sugar: A 5-membered ring with 3' hydroxyl, no 2' hydroxyl
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1O[C@@H]([C@H](O)C1)CO")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar backbone found (missing 5-membered ring with proper configuration)"
    
    # Nucleobase patterns for adenine, guanine, thymine, cytosine, uracil (varied configurations)
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2[nH]cnc2c1"),  # adenine-like
        Chem.MolFromSmarts("c1nc2[nH]c[nH]c2[nH]1"),  # guanine-like
        Chem.MolFromSmarts("c1c[nH]c(=O)nc1=O"),  # uracil-like
        Chem.MolFromSmarts("c1cc(=O)[nH]c[nH]1"),  # cytosine-like
        Chem.MolFromSmarts("c1c[nH]c(=O)nc1C")  # thymine-like adjustments
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns):
        return False, "No recognized nucleobase found (purine or pyrimidine derivatives)"

    # Check for a flexible 5'-phosphate group allowing different attachment types
    phosphate_pattern = Chem.MolFromSmarts("O[P](=O)(O)[O-]")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No flexible 5'-phosphate group match"

    return True, "Contains 2'-deoxyribose with 5'-phosphate and recognized nucleobase"