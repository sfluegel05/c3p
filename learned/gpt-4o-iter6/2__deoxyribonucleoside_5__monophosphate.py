"""
Classifies: CHEBI:18241 2'-deoxyribonucleoside 5'-monophosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

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

    # Check for deoxyribose sugar (without OH at 2' position)
    deoxyribose_pattern = Chem.MolFromSmarts("[C@@H]1O[C@H]([C@@H](COP(=O)([O-])O)[C@H]1O)[C@H]([CH2]O)O")
    if not mol.HasSubstructMatch(deoxyribose_pattern):
        return False, "No 2'-deoxyribose sugar structure found"

    # Check for phosphate group at 5' position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)([O-])O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No 5'-monophosphate group found"

    # Check for a nucleobase: either purine or pyrimidine pattern
    purine_pattern = Chem.MolFromSmarts("c1ncnc2n(cnc12)[C@H]")
    pyrimidine_pattern = Chem.MolFromSmarts("c1cc[nH]c(=O)n1")
    if not (mol.HasSubstructMatch(purine_pattern) or mol.HasSubstructMatch(pyrimidine_pattern)):
        return False, "No recognized purine or pyrimidine nucleobase found"

    return True, "Valid 2'-deoxyribonucleoside 5'-monophosphate with correct sugar, phosphate, and nucleobase"