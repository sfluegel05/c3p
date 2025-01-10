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

    # Define SMARTS for beta-D-glucuronic acid in its salt form, including stereochemistry variants
    glucuronate_smarts = [
        'O[C@H]1[C@@H](O)[C@H](O)[C@H](O[C@@H]1)C(=O)[O-]', # Major form
        'O[C@@H]1[C@@H](O)[C@H](O)[C@H](O[C@H]1)C(=O)[O-]', # Stereochemistry variants
        'O[C@H]1[C@H](O)[C@@H](O)[C@@H](O[C@@H]1)C(=O)[O-]'  # All possible isomers
    ]
    
    # Check for presence of the glucuronate moiety
    has_glucuronate = any(mol.HasSubstructMatch(Chem.MolFromSmarts(smarts)) for smarts in glucuronate_smarts)
    
    if not has_glucuronate:
        return False, "Beta-D-glucuronic acid moiety not found."

    # Comprehensive linkage patterns
    linkage_patterns = [
        '[O][C@H]',         # Ether link
        '[C](=O)[O][C]',    # Ester link
        '[N][C](=O)',       # Amide link
        '[O](S=O)',         # Sulfate link
        'c:[O][C]',         # Aromatic ether	      
        '[O][S](=O)(=O)[C]', # Sulfonation link
        '[C](=O)[N]',        # More generic amide
    ]

    # Check for any valid linkage pattern around the glucuronate moiety
    has_valid_linkage = any(mol.HasSubstructMatch(Chem.MolFromSmarts(link)) for link in linkage_patterns)
    
    if has_valid_linkage:
        return True, "Contains deprotonated beta-D-glucuronic acid moiety with valid linkage."
    
    return False, "Beta-D-glucuronic acid moiety found but no valid linkage detected."