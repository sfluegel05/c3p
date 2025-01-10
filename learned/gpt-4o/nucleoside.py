"""
Classifies: CHEBI:33838 nucleoside
"""
from rdkit import Chem

def is_nucleoside(smiles: str):
    """
    Determines if a molecule is a nucleoside based on its SMILES string.
    A nucleoside has both a nucleobase and either a ribose or deoxyribose sugar.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for more precise nucleobases
    adenine_pattern = Chem.MolFromSmarts('n1cnc2[nH]cnc2c1')  # Adenine
    guanine_pattern = Chem.MolFromSmarts('n1c2[nH]cnc2c(=O)[nH]c1=O')  # Guanine
    uracil_pattern = Chem.MolFromSmarts('n1c(=O)ccc[nH]1=O')  # Uracil
    thymine_pattern = Chem.MolFromSmarts('n1c(=O)ccn(C)c1=O')  # Thymine
    cytosine_pattern = Chem.MolFromSmarts('n1c(=O)ccn[nH]1')  # Cytosine

    nucleobase_patterns = [adenine_pattern, guanine_pattern, uracil_pattern, thymine_pattern, cytosine_pattern]

    # Check for nucleobases presence
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not has_nucleobase:
        return False, "No common nucleobase found"

    # Define more precise SMARTS patterns for ribose and deoxyribose sugar
    ribose_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H]([C@H](O)[C@H](O1)CO)O')  # Ribose
    deoxyribose_pattern = Chem.MolFromSmarts('O[C@H]1[C@@H]([C@H](O)[C@H](CO)O1)')  # Deoxyribose

    sugar_patterns = [ribose_pattern, deoxyribose_pattern]

    # Check for sugar presence
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"

    return True, "Contains both a nucleobase and either a ribose or deoxyribose sugar"