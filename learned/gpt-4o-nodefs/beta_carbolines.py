"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    
    A beta-carboline contains a pyrido[3,4-b]indole core with potential varied substitution.
    Enhance detection by looking for related core structures.
    
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
    
    # Define beta-carboline-like patterns potentially capturing variations
    beta_carboline_patterns = [
        Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3cccnc23"),  # 9H-pyrido[3,4-b]indole
        Chem.MolFromSmarts("c1ccc2c(c1)nc3[nH]ccc3c2"),  # Extended and substituted versions
        Chem.MolFromSmarts("c1ccc2c(c1)[nH]c3ccc[nH]c23"),  # N-substituted variants
        Chem.MolFromSmarts("c1ccc2c(c1)nc3cccnc23")  # Fully aromatic variation
    ]
    
    # Check if the molecule contains any of the beta-carboline core structures
    for pattern in beta_carboline_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains beta-carboline core structure or variant"
    
    return False, "Missing recognizable beta-carboline core structure"

# Example usage:
# result, reason = is_beta_carbolines("CN1C2=C(C=CC(=C2)OC)C3=C1[C@H](NCC34CCN(CC4)C(=O)C5=CC=C(C=C5)F)CO")
# print(result, reason)