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

    # Define a generalized SMARTS pattern for beta-carboline core (pyridoindole structure)
    # Allow for some variation in hydrogenation and substitutions
    pyridoindole_patterns = [
        Chem.MolFromSmarts("C12=CNC3=CC=CC=C3C(=N1)C=CC=C2"),  # Basic beta-carboline core
        Chem.MolFromSmarts("C1=C2CN=C(C=C2NC3=CC=CC=C31)"),     # Hydrogenated forms
        Chem.MolFromSmarts("n1cc2ccccn2c3ccccc13"),             # Fully aromatic variant
    ]

    # Check for matching the broad set of pyridoindole patterns
    for pattern in pyridoindole_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline core structure"

    return False, "Does not contain beta-carboline core structure"