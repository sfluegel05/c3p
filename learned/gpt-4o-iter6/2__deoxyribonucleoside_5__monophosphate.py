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

    # Patterns for 2'-deoxyribose, accounting for stereochemistry
    deoxyribose_patterns = [
        Chem.MolFromSmarts("O[C@H]1C[C@H](O)[C@@H](CO1)"),
        Chem.MolFromSmarts("O[C@@H]1C[C@H](O)[C@H](CO1)"),
        Chem.MolFromSmarts("O[C@H]1C[C@@H](O)[C@H](CO1)"),
        Chem.MolFromSmarts("O[C@@H]1C[C@@H](O)[C@H](CO1)")
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in deoxyribose_patterns):
        return False, "No valid 2'-deoxyribose sugar structure found"

    # Patterns for flexible 5'-monophosphate groups
    phosphate_patterns = [
        Chem.MolFromSmarts("COP(=O)(O)O"),
        Chem.MolFromSmarts("COP([O-])(=O)O"),
        Chem.MolFromSmarts("COP([O-])(=O)[O-]"),
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns):
        return False, "No valid 5'-monophosphate group found"

    # Patterns for comprehensive nucleobase representation, considering various types
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"),  # purine (e.g. adenine, guanine)
        Chem.MolFromSmarts("c1ccn(c2c(=O)n(cnc12)C)"),  # consistent pyrimidine
        Chem.MolFromSmarts("c1nc[nH]c(=O)n1"),  # further pyrimidine variations
        Chem.MolFromSmarts("c1ccc2[nH]c(nc2c1)C=O"),  # potential modifications
        Chem.MolFromSmarts("c1ncnc2[nH]cnc1c2")  # expanded purine configurations
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns):
        return False, "No recognized nucleobase found"

    return True, "Valid 2'-deoxyribonucleoside 5'-monophosphate with correct sugar, phosphate, and nucleobase"