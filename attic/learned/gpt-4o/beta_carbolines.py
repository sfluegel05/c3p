"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    A beta-carboline contains a pyridoindole core structure and can include 
    their hydrogenated derivatives.

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

    # Define the beta-carboline core pattern (pyridoindole structure)
    # Example SMARTS string may vary, it requires fine-tuning
    beta_carboline_pattern = Chem.MolFromSmarts('C1=NC2=CC=CC=C2C3=C1C=CC=C3')
    
    # Check for the beta-carboline pattern
    if not mol.HasSubstructMatch(beta_carboline_pattern):
        return False, "No beta-carboline structure found"

    return True, "Contains beta-carboline structure"

# Example usage:
smiles_example = "CC1=NC2=CC=CC=C2C3=C1C=CC=C3"  # Example beta-carboline SMILES
print(is_beta_carbolines(smiles_example))