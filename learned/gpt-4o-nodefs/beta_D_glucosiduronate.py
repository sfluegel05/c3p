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

    # SMARTS for beta-D-glucuronic acid moiety; ensure correct stereochemistry
    glucuronic_acid_pattern = Chem.MolFromSmarts("[C@@H]1(O[C@H]([C@@H]([C@H]([C@H]1O)O)O)C(=O)[O-])")
    if not mol.HasSubstructMatch(glucuronic_acid_pattern):
        return False, "Glucuronic acid moiety not found"

    # Check for ether or ester linkage patterns, allowing flexibility in attachment
    ether_linkage_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H]([C@H]([C@@H]([C@H]1O)O)O)C(=O)[O-]")
    ester_linkage_pattern = Chem.MolFromSmarts("O=C(O[C@H])1[C@H]([C@@H]([C@H]([C@H]1O)O)O)O")
    
    if not mol.HasSubstructMatch(ether_linkage_pattern) and not mol.HasSubstructMatch(ester_linkage_pattern):
        return False, "Glucuronic acid not attached via ether or ester linkage"

    return True, "Contains a beta-D-glucuronic acid moiety attached to a larger structure via ether or ester linkage"