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

    # Define a SMARTS pattern for ribose or deoxyribose
    sugar_pattern = Chem.MolFromSmarts("C1[C@H]([C@@H]([C@H](O1)COP(=O)(O)O)[OH])[OH]")
    
    # Define a SMARTS pattern for nucleobase (simplified examples)
    nucleobase_pattern = Chem.MolFromSmarts("n1cnc2c1ncn(c2)C=C")  # Adenine example

    # Define a SMARTS pattern for a phosphate group at the 5' position
    phosphate_pattern = Chem.MolFromSmarts("COP(=O)(O)O")
    
    # Check for the presence of ribose or deoxyribose
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No ribose or deoxyribose sugar found"

    # Check for the presence of a nucleobase
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"

    # Check for the presence of a phosphate group at the 5' position
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "No phosphate group found at the 5' position"

    return True, "Contains sugar, nucleobase, and phosphate group at the 5' position"