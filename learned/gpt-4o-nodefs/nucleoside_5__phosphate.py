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

    # Flexible patterns for ribose or deoxyribose sugars with phosphate
    sugar_phosphate_pattern = Chem.MolFromSmarts("C1C(O[C@H]2CO[C@@H](C2)O)C(COP(O)(=O)[O])O1")

    # Flexible SMARTS patterns for canonical nucleobases
    adenine_pattern = Chem.MolFromSmarts("n1cnc2ncnc2n1")
    guanine_pattern = Chem.MolFromSmarts("n1c2c(ncnc2n(c1))O")
    cytosine_pattern = Chem.MolFromSmarts("n1c2c(ncnc2n(c1))N")
    thymine_pattern = Chem.MolFromSmarts("C1=CN(C=O)C(=O)N(C1)=O")
    uracil_pattern = Chem.MolFromSmarts("C1=CN(C=O)C(=O)NC1")

    nucleobase_patterns = [adenine_pattern, guanine_pattern, cytosine_pattern, thymine_pattern, uracil_pattern]

    # Check for the presence of a ribose or deoxyribose sugar with phosphate
    if not mol.HasSubstructMatch(sugar_phosphate_pattern):
        return False, "No 5'-phosphate sugar found"

    # Check for the presence of any nucleobase
    if not any(mol.HasSubstructMatch(base) for base in nucleobase_patterns):
        return False, "No nucleobase found"

    return True, "Contains sugar, nucleobase, and phosphate group at the 5' position"