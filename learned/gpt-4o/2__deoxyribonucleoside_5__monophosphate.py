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
    
    # Pattern for 2'-deoxyribose sugar with no hydroxyl at 2' position, hydroxy at 3', and correct stereochemistry
    deoxyribose_pattern = Chem.MolFromSmarts("C1([C@@H](O)C([C@H](O)CO[P](=O)(O)O)O[C@H]1)")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar backbone found"

    # Nucleobase patterns for adenine, guanine, thymine, cytosine, and uracil
    purine_base_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2")
    pyrimidine_base_pattern = Chem.MolFromSmarts("n1ccn(C)c(=O)c1")
    
    is_purine = mol.HasSubstructMatch(purine_base_pattern)
    is_pyrimidine = mol.HasSubstructMatch(pyrimidine_base_pattern)
    
    if not (is_purine or is_pyrimidine):
        return False, "No nucleobase found"

    # Check for a phosphate group at 5' location
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No 5'-phosphate group found"

    return True, "Contains 2'-deoxyribose with 5'-phosphate and nucleobase"