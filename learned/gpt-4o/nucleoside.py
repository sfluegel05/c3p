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
        bool: True if the molecule is a nucleoside, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for nucleobases (include variations and possible tautomers)
    adenine_pattern = Chem.MolFromSmarts('n1c[nH]c2c1ncnc2')  # Adenine
    guanine_pattern = Chem.MolFromSmarts('c1c([nH]c2nc[nH]c2n1)N')  # Guanine
    uracil_pattern = Chem.MolFromSmarts('c1cnc[nH]c1=O')  # Uracil
    thymine_pattern = Chem.MolFromSmarts('c1cnc[nH]c1=O')  # Thymine
    cytosine_pattern = Chem.MolFromSmarts('c1cncnc1=O')  # Cytosine

    nucleobase_patterns = [adenine_pattern, guanine_pattern, uracil_pattern, thymine_pattern, cytosine_pattern]

    # Check for presence of a nucleobase
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in nucleobase_patterns)
    if not has_nucleobase:
        return False, "No common nucleobase found"

    # Define SMARTS patterns for ribose and deoxyribose sugar
    ribose_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H]([C@H](O)[C@@H](O1)CO)')  # Ribose
    deoxyribose_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H]([C@@H](O)[C@H](CO)O1)')  # Deoxyribose

    sugar_patterns = [ribose_pattern, deoxyribose_pattern]

    # Check for presence of ribose or deoxyribose sugar
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugar_patterns)
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"

    return True, "Contains both a nucleobase and either a ribose or deoxyribose sugar"