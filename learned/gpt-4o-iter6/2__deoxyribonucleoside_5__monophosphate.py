"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    This class involves a 2'-deoxyribose sugar with a 5'-phosphate and one of the typical nucleobases.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2'-deoxyribonucleoside 5'-monophosphate, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More precise 2'-deoxyribose pattern considering stereochemistry
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]1([C@H]([C@@H]([C@@H](O1)CO)O)O)")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar structure found"
    
    # Flexible pattern for monophosphate group, allowing for charge variations
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        phosphate_pattern_charged = Chem.MolFromSmarts("COP([O-])(=O)O")
        if not mol.HasSubstructMatch(phosphate_pattern_charged):
            return False, "No 5'-monophosphate group found"
    
    # Comprehensive nucleobase patterns
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),  # purines (generic)
        Chem.MolFromSmarts("n1cnc2ncnc(c12)N"),  # adenine-like
        Chem.MolFromSmarts("n1c[nH]cnc1")  # guanine-like
    ] + [
        Chem.MolFromSmarts("c1cnc(=O)nc1"),  # pyrimidine
        Chem.MolFromSmarts("n1cc(cn1C)C=O"), # thymine-like
        Chem.MolFromSmarts("c1[nH]cnc1=N"),  # cytosine-like
        Chem.MolFromSmarts("[nH]1c(=O)[nH]c(=O)c([nH]1)N")  # uracil-like
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns):
        return False, "No recognized purine or pyrimidine nucleobase found"
    
    return True, "Valid 2'-deoxyribonucleoside 5'-monophosphate with correct sugar, phosphate, and nucleobase"