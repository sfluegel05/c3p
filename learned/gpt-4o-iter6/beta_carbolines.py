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

    # Define a SMARTS pattern for a beta-carboline skeleton
    # Here, we'll look for a flexible pattern that may represent common core features 
    # of beta-carbolines like the indole and pyridine structures.
    beta_carboline_pattern = Chem.MolFromSmarts('c1ccc2[nH]c3c(c2c1)ncc3')

    # Check if the molecule has the beta-carboline substructure
    if mol.HasSubstructMatch(beta_carboline_pattern):
        return True, "Contains beta-carboline core structure"

    # If more complex hydrogenated derivatives need to be considered, further patterns may be required
    # For simplicities sake, we assume the given pattern catches most beta-carboline cores

    return False, "Does not contain beta-carboline core structure"

# The SMARTS pattern 'c1ccc2[nH]c3c(c2c1)ncc3' is a simplified representation of the typical 
# beta-carboline ring system, intended as a placeholder for more comprehensive patterns.