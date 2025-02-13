"""
Classifies: CHEBI:83411 beta-D-glucosiduronate
"""
"""
Classifies: CHEBI:10004 beta-D-glucosiduronate
"""
from rdkit import Chem

def is_beta_D_glucosiduronate(smiles: str):
    """
    Determines if a molecule is a beta-D-glucosiduronate based on its SMILES string.
    A beta-D-glucosiduronate contains a deprotonated β-D-glucuronic acid moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for beta-D-glucuronic acid moiety in deprotonated form
    glucuronate_pattern = Chem.MolFromSmarts('[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(=O)[O-]')
    if not mol.HasSubstructMatch(glucuronate_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    # Check attachment type (e.g., it is often via an ether linkage)
    # This could be broadened based on additional requirements or examples
    ether_linkage_pattern = Chem.MolFromSmarts('O[C@@H]1O[C@@H]([C@@H](O)[C@H](O)[C@H]1O)C(=O)[O-]')
    if not mol.HasSubstructMatch(ether_linkage_pattern):
        return False, "Beta-D-glucuronic acid moiety not properly attached"

    return True, "Contains deprotonated beta-D-glucuronic acid moiety in correct orientation"