"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains a sugar moiety (often a pyranose or furanose ring or open structure)
    with one or more phosphate groups attached.

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

    # Flexible wildcard patterns for phosphate groups, accommodating different states and linkages
    phosphate_patterns = [
        Chem.MolFromSmarts("[OP](=O)(O)O"),     # General protonated phosphate
        Chem.MolFromSmarts("[OP](=O)(O)[O-]"), # General deprotonated phosphate
        Chem.MolFromSmarts("[OP](=O)(O)C"),    # Phosphate ester linkage
        Chem.MolFromSmarts("[OP](=O)(O)S")     # Phosphorothioate
    ]

    # Revised and more comprehensive sugar ring and carbohydrates patterns
    sugar_patterns = [
        Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1"),               # Pyranose form
        Chem.MolFromSmarts("C1OC(O)C(O)C1"),                    # Furanose form
        Chem.MolFromSmarts("O[C@@H]1[C@@H](O)[C@H](O)C(O)C1"),  # Ribose/Aldose
        Chem.MolFromSmarts("[C@H](O)[C@@H](O)C(O)(O)C=O"),     # Open chain sugars (like glucose, fructose)
    ]

    # Check for at least one phosphate group pattern match
    has_phosphate = any(mol.HasSubstructMatch(pat) for pat in phosphate_patterns)
    if not has_phosphate:
        return False, "No phosphate group found"

    # Revised broader checking for presence of any sugar-related structure
    has_sugar_structure = any(mol.HasSubstructMatch(pat) for pat in sugar_patterns)
    if not has_sugar_structure:
        return False, "No sugar structure detected"

    # If both sugar and phosphate are present, classify as phospho sugar
    return True, "Contains a phosphate group attached to a sugar moiety"