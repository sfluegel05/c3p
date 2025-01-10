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

    # 2'-deoxyribose sugar pattern with consideration to stereochemistry
    deoxyribose_pattern = Chem.MolFromSmarts("O[C@H]1[C@H]([C@@H](O)[C@@H](C1)O)")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar structure found"
    
    # Flexible pattern for monophosphate group with variations
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    phosphate_pattern_charged = Chem.MolFromSmarts("COP([O-])(=O)O")
    if not mol.HasSubstructMatch(phosphate_pattern) and not mol.HasSubstructMatch(phosphate_pattern_charged):
        return False, "No 5'-monophosphate group found"

    # Comprehensive nucleobase patterns, including alternative configurations
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),  # generic purine
        Chem.MolFromSmarts("c1ncnc(=N1)N"),  # adenine-like alternative
        Chem.MolFromSmarts("n1c(nc2n(c1)nc2N)"),  # guanine-like alternative

        Chem.MolFromSmarts("c1ncnc(=O)n1"),  # generic pyrimidine
        Chem.MolFromSmarts("n1c(cnc1O)C=O"), # thymine-like alternative
        Chem.MolFromSmarts("n1cc(cn1)[NH]"),  # cytosine-like alternative
        Chem.MolFromSmarts("n1c(=O)c(=O)nc[nH]1"),  # uracil-like alternative
    ]

    if not any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns):
        return False, "No recognized purine or pyrimidine nucleobase found"

    return True, "Valid 2'-deoxyribonucleoside 5'-monophosphate with correct sugar, phosphate, and nucleobase"