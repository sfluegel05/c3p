"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains a sugar moiety (often as a furanose or pyranose)
    with one or more attached phosphate groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phospho sugar, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a phosphate group pattern (including deprotonated variation)
    phosphate_patterns = [
        Chem.MolFromSmarts("OP(=O)(O)O"),  # Regular phosphate
        Chem.MolFromSmarts("OP(=O)([O-])[O-]"),  # Deprotonated phosphates
    ]
    
    # Check for at least one phosphate group pattern match
    has_phosphate = any(mol.HasSubstructMatch(pat) for pat in phosphate_patterns)
    if not has_phosphate:
        return False, "No phosphate group found"

    # Define pyranose and furanose ring patterns (common sugar structures)
    pyranose_pattern = Chem.MolFromSmarts("O1CCCC1[OH]")
    furanose_pattern = Chem.MolFromSmarts("O1CCC1[OH]")

    # Check for sugar structure (pyranose or furanose)
    has_sugar_ring = mol.HasSubstructMatch(pyranose_pattern) or mol.HasSubstructMatch(furanose_pattern)
    if not has_sugar_ring:
        return False, "No sugar ring structure detected"

    # If both sugar and phosphate are present, classify as phospho sugar
    return True, "Contains a phosphate group attached to a sugar moiety (furanose or pyranose ring)"