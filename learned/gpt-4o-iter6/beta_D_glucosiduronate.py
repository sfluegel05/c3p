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

    # Improved SMARTS pattern for beta-D-glucuronic acid moiety in deprotonated form
    # Consider common patterns for glucuronic acid moiety with oxygen and carbon orientation
    glucuronate_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H](O)[C@@H](O)[C@H](O[C@H]1)C(=O)[O-]')
    if not mol.HasSubstructMatch(glucuronate_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    # Add patterns to detect ether or ester linkage in varying scenarios
    possible_linkages = [
        'O[C@]',  # common ether link
        'O=C-O[#6]',  # ester linkage possibly present between glucuronide and aglycone
    ]
    for linkage in possible_linkages:
        linkage_pattern = Chem.MolFromSmarts(linkage)
        if mol.HasSubstructMatch(linkage_pattern):
            return True, "Contains deprotonated beta-D-glucuronic acid moiety with correct linkage"

    return False, "Beta-D-glucuronic acid moiety not properly attached"