"""
Classifies: CHEBI:16701 nucleoside 5'-phosphate
"""
from rdkit import Chem

def is_nucleoside_5__phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside 5'-phosphate based on its SMILES string.
    A nucleoside 5'-phosphate consists of a sugar, a nucleobase, and a phosphate group at the 5' position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a nucleoside 5'-phosphate, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for ribose or deoxyribose
    sugar_pattern = Chem.MolFromSmarts("C1OC([C@@H](O)C(O)C1)COP(O)(=O)[O]")  # Can match both ribose and deoxyribose
    
    # Define SMARTS patterns for canonical nucleobases
    adenine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")
    guanine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2O")
    cytosine_pattern = Chem.MolFromSmarts("n1cnc2c1ncnc2N")
    thymine_pattern = Chem.MolFromSmarts("C1=CN(C=O)C(=O)NC1=O")
    uracil_pattern = Chem.MolFromSmarts("C1=CN(C=O)C(=O)NC1=O")

    nucleobase_patterns = [adenine_pattern, guanine_pattern, cytosine_pattern, thymine_pattern, uracil_pattern]

    # Define a SMARTS pattern for a phosphate group at the 5' position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    
    # Check for the presence of ribose or deoxyribose
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose or deoxyribose sugar found"

    # Check for the presence of any nucleobase
    if not any(mol.HasSubstructMatch(base) for base in nucleobase_patterns):
        return False, "No nucleobase found"

    # Check for the presence of a phosphate group at the 5' position
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found at the 5' position"

    return True, "Contains sugar, nucleobase, and phosphate group at the 5' position"