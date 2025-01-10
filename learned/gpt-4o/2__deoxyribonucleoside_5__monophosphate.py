"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.
    A 2'-deoxyribonucleoside monophosphate typically contains a deoxyribose sugar, 
    a nitrogenous base, and a phosphate group at the 5' position.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule fits the class, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for 2'-deoxyribose sugar pattern
    # Revising pattern to consider only the structural part and not stereochemistry.
    # 5-membered ring with oxygen (furanose)
    deoxyribose_pattern = Chem.MolFromSmarts("C1C(C(C(C1O)O)O)CO")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose structure found"

    # Check for phosphate group pattern at 5' position
    # Broadening pattern to match potential variations including charged forms
    phosphate_patterns = [
        Chem.MolFromSmarts("COP(=O)(O)[O-]"),  # charged 
        Chem.MolFromSmarts("COP(=O)(O)O"),    # neutral phosphate
        Chem.MolFromSmarts("COP(=O)(O)OC"),   # ester link
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in phosphate_patterns):
        return False, "No 5'-phosphate group found"

    # Check for typical nucleobases pattern - covering more permutations
    base_patterns = [
        Chem.MolFromSmarts("n1cnc2c1ncnc2"), # adenine or purine bases
        Chem.MolFromSmarts("n1c(O)nc2c1ncnc2"), # oxo-purines
        Chem.MolFromSmarts("n1ccnc1=O"), # pyrimidine bases
    ]
    if not any(mol.HasSubstructMatch(base) for base in base_patterns):
        return False, "No typical nucleobase structure found"

    return True, "Structure fits 2'-deoxyribonucleoside 5'-monophosphate definition"