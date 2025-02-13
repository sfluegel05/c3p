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

    # Define a set of patterns to capture the beta-carboline and its derivatives
    # This includes basic and hydrogenated versions
    beta_carboline_patterns = [
        Chem.MolFromSmarts("c1cc2c[nH]c3cccc(c3n2)c1"),    # Standard beta-carboline core
        Chem.MolFromSmarts("c1c2[nH]c3c([nH]2)cccc3c1"),   # Hydrogenated variant
        Chem.MolFromSmarts("C1Cc2c([nH]c3c2cccc3)NC=C1")   # Another common hydrogenation variant
    ]
    
    # Check if the molecule matches any of the beta-carboline patterns
    for pattern in beta_carboline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline skeleton or derivative"

    return False, "Does not contain beta-carboline skeleton"