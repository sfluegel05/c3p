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
    # Adjusted to match the primary sugar pattern with C6 deprotonated
    glucuronate_pattern = Chem.MolFromSmarts('O[C@@H]1[C@H](O)[C@H](O)[C@@H](O[C@H]1[C@@H]=O)C(=O)[O-]')
    if not mol.HasSubstructMatch(glucuronate_pattern):
        return False, "No beta-D-glucuronic acid moiety found"

    # Broadened pattern acceptance for linkage - consider various linkages
    possible_linkages = [
        'O[C@]',  # Ether
        # These can be expanded if other types of linkages are identified as characteristic
    ]
    for linkage in possible_linkages:
        linkage_pattern = Chem.MolFromSmarts(linkage)
        if mol.HasSubstructMatch(linkage_pattern):
            return True, "Contains deprotonated beta-D-glucuronic acid moiety in correct orientation"

    return False, "Beta-D-glucuronic acid moiety not properly attached"