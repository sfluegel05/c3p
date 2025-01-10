"""
Classifies: CHEBI:33447 phospho sugar
"""
from rdkit import Chem

def is_phospho_sugar(smiles: str):
    """
    Determines if a molecule is a phospho sugar based on its SMILES string.
    A phospho sugar contains a sugar moiety (often a pyranose or furanose ring or open chain)
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

    # Expand patterns to include more cases and states of phosphate groups
    phosphate_patterns = [
        Chem.MolFromSmarts("P(=O)(O)O"),      # General phosphate
        Chem.MolFromSmarts("P(=O)(O)[O-]"),   # Deprotonated phosphate
        Chem.MolFromSmarts("OP(=O)(O)O"),     # Phosphate ester
        Chem.MolFromSmarts("OP(=O)([O-])O")   # Phosphate in anionic form
    ]

    # Expand sugar detection to cover open chain and various configurations
    sugar_patterns = [
        Chem.MolFromSmarts("OC[C@H]1O[C@H](O)[C@@H](O)[C@H](O)[C@H]1O"),  # Pyranose alpha
        Chem.MolFromSmiles("C1=COC(O)C(O)C1O"),                          # Furanose
        Chem.MolFromSmiles("OC[C@H](O)[C@H](O)[C@H](O)C=O"),              # Aldose open chain
        Chem.MolFromSmarts("C(C(C(C=O)O)O)O")                            # Generic open form
    ]

    # Checking for at least one phosphate group pattern match
    has_phosphate = any(mol.HasSubstructMatch(pat) for pat in phosphate_patterns)
    if not has_phosphate:
        return False, "No phosphate group found"

    # Checking for presence of sugar-related structure
    has_sugar_structure = any(mol.HasSubstructMatch(pat) for pat in sugar_patterns)
    if not has_sugar_structure:
        return False, "No sugar structure detected"

    # If both sugar and phosphate are present, classify as phospho sugar
    return True, "Contains a phosphate group attached to a sugar moiety"