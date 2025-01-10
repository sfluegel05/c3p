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
    
    # SMARTS pattern for ribose or deoxyribose with variable phosphate attachment at 5'
    sugar_phosphate_pattern = Chem.MolFromSmarts("C1C(O)C(O)C(O)C1O[*:1]P(=O)(O)O")
    
    # Nucleobase SMARTS patterns: expanded to consider variations (the generic form of pyrimidine, purine, including tautomeric shifts)
    nucleobase_patterns = [
        Chem.MolFromSmarts("n1cnc2ncnc2n1"),  # Adenine akin, including modifications
        Chem.MolFromSmarts("n1c2c(ncnc2n(c1))O"),  # Guanine akin
        Chem.MolFromSmarts("n1cnc2[nH]cnc2n1"),  # Cytosine, working-bearing nitrogenous atoms
        Chem.MolFromSmarts("C1=C[NH]C(=O)N(C1)=O"),  # Thymine/uracil reflection
        # More could be added for less common modifications.
    ]

    # Check for the presence of a sugar with phosphate group
    if not mol.HasSubstructMatch(sugar_phosphate_pattern):
        return False, "No valid 5'-phosphate sugar recognized"
    
    # Check for the presence of a nucleobase
    if not any(mol.HasSubstructMatch(base) for base in nucleobase_patterns):
        return False, "No nucleobase found"

    return True, "Contains a nucleobase and 5'-phosphate sugar."