"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined SMARTS for beta-D-glucuronic acid moiety with exact stereochemistry
    glucuronic_acid_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "Glucuronic acid moiety not correctly found or stereochemistry missed"

    # Comprehensive check for ether and ester linkages incorporating flexibility
    # More patterns to ensure coverage of varied connections
    # Note: Using patterns that ensure the glucuronic acid is linked appropriately
    ether_linkage_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]~*")  # Allowing any atom after O-
    ester_linkage_pattern = Chem.MolFromSmarts("O=C(O[C@H]1[C@H]([C@@H]([C@H]([C@H]1O)O)O)O)~*")  # Again, allowing connections

    # Both or either should match as per the linkage requirements; adding more flexibility
    if not (mol.HasSubstructMatch(ether_linkage_pattern) or mol.HasSubstructMatch(ester_linkage_pattern)):
        return False, "Glucuronic acid linkage not identified; revise ether/ester condition matching"

    return True, "Contains a beta-D-glucuronic acid moiety linked via recognized ether or ester connection"