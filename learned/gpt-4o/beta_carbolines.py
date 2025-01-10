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
    
    # Define a more flexible beta-carboline pattern
    # The pattern should match a pyridoindole core, with allowance for partial unsaturation
    beta_carboline_pattern = Chem.MolFromSmarts('n1c2ccc(cc2c3c1cccc3)')  # More flexible pyridoindole

    # Check for the beta-carboline pattern
    if not mol.HasSubstructMatch(beta_carboline_pattern):
        return False, "No beta-carboline structure found"

    return True, "Contains a beta-carboline structure"

# Example usage with one of the provided example SMILES
smiles_example = "CCCNC(=O)N1CC2(C1)CN([C@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)S(=O)(=O)C"
print(is_beta_carbolines(smiles_example))