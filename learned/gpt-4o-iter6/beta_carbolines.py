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

    # Define enhanced SMARTS pattern for beta-carboline core
    # capturing more structural variations found in example structures
    beta_carboline_patterns = [
        Chem.MolFromSmarts("n1c2c(cccc2)c3c1nccc3"),  # Basic beta-carboline core
        Chem.MolFromSmarts("n1c2c(cccc2)c3c1nc[nH]c3"),  # Various nitrogen placements
        Chem.MolFromSmarts("n1c2c(ccncc2)c3c1nccc3"),  # More nitrogen substitutions
        Chem.MolFromSmarts("n1c2c(cccnc2)c3c1nc[nH]c3"),  # Dihydro variants
        Chem.MolFromSmarts("C1Cc2c(C1)nc3c(n2)cccc3"),  # Partial hydrogenation
        Chem.MolFromSmarts("c1c[nH]c2c1ncc3c2cccn3"),  # Indole-like structures
    ]

    # Check for matching any of the expanded beta-carboline patterns
    for pattern in beta_carboline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline core structure"

    return False, "Does not contain beta-carboline core structure"