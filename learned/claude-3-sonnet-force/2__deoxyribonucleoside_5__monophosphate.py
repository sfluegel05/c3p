"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
"""
Classifies: CHEBI:18115 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_2__deoxyribonucleoside_5__monophosphate(smiles: str):
    """
    Determines if a molecule is a 2'-deoxyribonucleoside 5'-monophosphate based on its SMILES string.

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

    # Define SMARTS patterns for components
    deoxyribose_pattern = Chem.MolFromSmarts("[C@H]([C@@H]([C@H]([C@H](CO)O)O)O)[C@@H]1[CH]O1")
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(=O)O")
    adenine_pattern = Chem.MolFromSmarts("nc1ncnc2n(cnc12)N")
    guanine_pattern = Chem.MolFromSmarts("Nc1nc2nc([nH]c2c(=O)[nH]1)N")
    cytosine_pattern = Chem.MolFromSmarts("Nc1cc[nH]c(=O)n1")
    thymine_pattern = Chem.MolFromSmarts("Cc1cn([C@H]2C[C@H](O)[C@@H](COP(O)(O)=O)O2)c(=O)[nH]c1=O")
    uracil_pattern = Chem.MolFromSmarts("O=c1cc[nH]c(=O)[nH]1")

    # Check for 2'-deoxyribose sugar
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar found"

    # Check for phosphate group
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found"

    # Check for nucleobase
    nucleobase_found = False
    for pattern in [adenine_pattern, guanine_pattern, cytosine_pattern, thymine_pattern, uracil_pattern]:
        if mol.HasSubstructMatch(pattern):
            nucleobase_found = True
            break
    if not nucleobase_found:
        return False, "No valid nucleobase found"

    # Check for connection between deoxyribose and phosphate
    deoxyribose_phosphate_pattern = Chem.MolFromSmarts("[C@H]([C@@H]([C@H]([C@H](CO)O)O)O)[C@@H]1[CH]O1OP(O)(=O)O")
    if not mol.HasSubstructMatch(deoxyribose_phosphate_pattern):
        return False, "Deoxyribose and phosphate not connected"

    # Check for connection between deoxyribose and nucleobase
    deoxyribose_nucleobase_pattern = Chem.MolFromSmarts("[C@H]([C@@H]([C@H]([C@H](CO)O)O)O)[C@@H]1[CH]O1~*")
    if not mol.HasSubstructMatch(deoxyribose_nucleobase_pattern):
        return False, "Deoxyribose and nucleobase not connected"

    return True, "Contains 2'-deoxyribose sugar with 5'-phosphate group and a valid nucleobase"