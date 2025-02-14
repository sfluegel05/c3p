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
    # Patterns include a broader range for hydrogenated derivatives and possible ring closures.
    beta_carboline_patterns = [
        Chem.MolFromSmarts("c1cc2[nH]cc(c3ccccc32)c1"),  # Standard beta-carboline core
        Chem.MolFromSmarts("c1cnc2c3ccccc3[nH]c2c1"),    # Classic and some hydrogenated forms
        Chem.MolFromSmarts("c1[nH]cc2c3ccccc3cnc2c1"),   # Another common hydrogenation
        Chem.MolFromSmarts("C1=CC=2C([nH]c3ccccc23)=NC1") # Variations of pyridoindole with shifts
    ]

    # Check if the molecule matches any of the beta-carboline patterns
    for pattern in beta_carboline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline skeleton or derivative"

    return False, "Does not contain beta-carboline skeleton"