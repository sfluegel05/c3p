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
    A beta-D-glucosiduronate contains a deprotonated beta-D-glucuronic acid moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-D-glucosiduronate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."

    # Define more comprehensive SMARTS patterns for beta-D-glucuronic acid, considering various stereochemistry configurations
    glucuronate_smarts = [
        'O[C@H]1[C@@H](O)[C@H](O)[C@H](O[C@@H]1)C(=O)[O-]',  # Normal anionic form
        'O[C@@H]1[C@@H](O)[C@H](O)[C@H](O[C@H]1)C(=O)[O-]', # Reverse chiral centers
        'O[C@H]1[C@H](O)[C@@H](O)[C@@H](O[C@@H]1)C(=O)[O-]'  # Any chiral center possibilities
    ]
    
    # Check for the presence of beta-D-glucuronic acid moiety using the SMARTS patterns
    has_glucuronate = any(mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for smarts in glucuronate_smarts)
    
    if not has_glucuronate:
        return False, "No beta-D-glucuronic acid moiety found."

    # Define potential linkage patterns for linkage types around the glucuronate moiety
    linkage_patterns = [
        '[O][#6]',         # Any carbon to oxygen linkage, possible ester or ether
        '[C](=O)[O][C]',   # Ester link pattern with carbonyl group
        '[N][C](=O)',      # Amide linkage pattern
        '[O][S](=O)(=O)',  # Sulfate linkage consideration
        'c:[O][C]',        # Aromatic ethers
        '[C](=O)[N]',      # Possible generic amidic linkage
        '[O]S(C)[C]',      # Specifically looking into sulfation type linkages
    ]

    # Check the molecule for any of the valid linkage types
    has_valid_linkage = any(mol.HasSubstructMatch(Chem.MolFromSmarts(link)) for link in linkage_patterns)
    
    if has_valid_linkage:
        return True, "Contains deprotonated beta-D-glucuronic acid moiety with valid linkage."
    
    return False, "Beta-D-glucuronic acid moiety found but lacks recognized valid linkages."