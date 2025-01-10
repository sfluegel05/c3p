"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline is defined as any pyridoindole containing a beta-carboline skeleton 
    and their hydrogenated derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define enhanced SMARTS patterns for beta-carboline core 
    # to accommodate different levels of saturation and substitutions
    beta_carboline_patterns = [
        Chem.MolFromSmarts("n1c2c(c[nH]c2)c3c1cccc3"),  # Essential beta-carboline core
        Chem.MolFromSmarts("n1c2cccc3c2c(cn1)c[nH]3"),  # Alternative arrangement
        Chem.MolFromSmarts("n1c2ccccc2c3c1c[nH]cn3"),  # More complex core with saturation
        Chem.MolFromSmarts("n1c2cc[nH]c2c3c1cccn3"),   # Including partially hydrogenated structures
    ]

    # Check for matching any of the expanded beta-carboline patterns
    for pattern in beta_carboline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline core structure"
    
    return False, "Does not contain beta-carboline core structure"