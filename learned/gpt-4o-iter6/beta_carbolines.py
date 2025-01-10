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
        return None, "Invalid SMILES string"

    # Define SMARTS patterns for beta-carboline skeleton and variants
    # Basic beta-carboline core (pyridoindole structure)
    pyridoindole_pattern = Chem.MolFromSmarts('C12=CNC3=CC=CC=C3C(=N1)C=CC=C2')  # Hydride variants

    # Check for basic pyridoindole patterns and their hydrogenated derivatives
    if mol.HasSubstructMatch(pyridoindole_pattern):
        return True, "Contains beta-carboline core structure"

    # Additional specific checks for other structural variants could be added 
    # based on more detailed examination of known beta-carbolines.

    return False, "Does not contain beta-carboline core structure"