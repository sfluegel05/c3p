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
        return (None, "Invalid SMILES string")

    # Define SMARTS patterns for nucleobases
    adenine_pattern = Chem.MolFromSmarts('c1ncnc2n(cnc12)')
    guanine_pattern = Chem.MolFromSmarts('c1nc2c(nc[nH]c2n1)')
    uracil_pattern = Chem.MolFromSmarts('c1cc(=O)[nH]c(=O)n1')
    thymine_pattern = Chem.MolFromSmarts('c1c[nH]c(=O)n1C')
    cytosine_pattern = Chem.MolFromSmarts('c1c[nH]c(=O)nc1')
    
    nucleobases_patterns = [adenine_pattern, guanine_pattern,
                            uracil_pattern, thymine_pattern, cytosine_pattern]
    
    # Check for nucleobases presence
    has_nucleobase = any(mol.HasSubstructMatch(pattern) for pattern in nucleobases_patterns)
    if not has_nucleobase:
        return False, "No common nucleobase found"

    # Define SMARTS pattern for a ribose or deoxyribose sugar
    ribose_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](O)[C@H](CO)O1')  # Ribofuranose
    deoxyribose_pattern = Chem.MolFromSmarts('O[C@H]1[C@H](O)[C@@H](CO)C[C@H]1O')  # Deoxyribofuranose
    
    sugars_patterns = [ribose_pattern, deoxyribose_pattern]
    
    # Check for sugar presence
    has_sugar = any(mol.HasSubstructMatch(pattern) for pattern in sugars_patterns)
    if not has_sugar:
        return False, "No ribose or deoxyribose sugar found"

    return True, "Contains both a nucleobase and either a ribose or deoxyribose sugar"