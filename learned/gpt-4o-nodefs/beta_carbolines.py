"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.

    A beta-carboline contains a pyrido[3,4-b]indole core with potential varied substitution.
    Enhanced detection by looking for related core structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a beta-carboline, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES to rdkit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define beta-carboline-like patterns that potentially capture variations
    # Using SMARTS to match a broader variety of beta-carboline cores
    beta_carboline_patterns = [
        Chem.MolFromSmarts("c1nc2ccc3c(c2c(c1)[nH]3)"),  # Basic 9H-pyrido[3,4-b]indole structure
        Chem.MolFromSmarts("c1nc2cccc3[nH]c(c1)c2cc3"),  # Variants with changes in the aromatic ring
        Chem.MolFromSmarts("c1ccc2nc3[nH]cc(c3c2c1)"),  # Tautomeric and alternative resonance forms
        Chem.MolFromSmarts("c1cc2ccc3[nH]c(c2c(c1)n3)"),  # Extended 5-membered ring systems
        Chem.MolFromSmarts("c1ccc2nc3ccc(c3[nH]2)c1"),  # Checking for various substitutions at nitrogen
    ]
    
    # Check if the molecule contains any of the beta-carboline core structures
    for pattern in beta_carboline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline core structure or variant"
    
    return False, "Missing recognizable beta-carboline core structure"

# Example usage:
# result, reason = is_beta_carbolines("CN1C2=C(C=CC(=C2)OC)C3=C1[C@H](NCC34CCN(CC4)C(=O)C5=CC=C(C=C5)F)CO")
# print(result, reason)