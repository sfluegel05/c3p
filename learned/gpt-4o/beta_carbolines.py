"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is defined as a pyridoindole containing a beta-carboline skeleton
    and their hydrogenated derivatives.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for the beta-carboline skeleton
    # These patterns must account for various configurations, and hydrogenated derivatives
    beta_carboline_patterns = [
        Chem.MolFromSmarts("c1cc2c3ccccc3nc2[nH]1"),   # Classic beta-carboline core
        Chem.MolFromSmarts("c1ncc2c3ccccc3[nH]c2c1"),  # Hydrogenated variant
        Chem.MolFromSmarts("C1=CC2=C3C(=CC=CC3=NC2=CN1)[nH]"),  # Another common structure
    ]

    # Check if the molecule matches any of the beta-carboline patterns
    for pattern in beta_carboline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline skeleton or derivative"

    return False, "Does not contain beta-carboline skeleton"