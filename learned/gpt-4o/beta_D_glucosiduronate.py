"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate features a beta-D-glucuronic acid moiety deprotonated
    at the carboxyl group, typically connected via an O-linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-D-glucosiduronate; otherwise False
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # More flexible SMARTS pattern for beta-D-glucuronic acid moiety
    beta_d_glucuronic_acid_pattern = Chem.MolFromSmarts("O[C@H]1[C@@H](O[C@H](O)[C@@H](O)[C@H]1)C(=O)[O-]")

    if not mol.HasSubstructMatch(beta_d_glucuronic_acid_pattern):
        return False, "No beta-D-glucuronic acid moiety found"
    
    # Check for O-linkage pattern ensuring sugar moiety contributes in a linkage
    o_linkage_pattern = Chem.MolFromSmarts("[C@H]1[C@@H](O[C@H](O)[C@@H](O)[C@H]1)C(=O)[O-]O[*]")
    
    if not mol.HasSubstructMatch(o_linkage_pattern):
        return False, "No appropriate O-linkage found with glucuronic acid moiety"

    return True, "Contains beta-D-glucuronic acid moiety with O-linkage and proper deprotonation"