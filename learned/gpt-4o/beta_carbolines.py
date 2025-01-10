"""
Classifies: CHEBI:60834 beta-carbolines
"""
from rdkit import Chem

def is_beta_carbolines(smiles: str):
    """
    Determines if a molecule is a beta-carboline based on its SMILES string.
    Beta-carbolines have a pyridoindole core and can have hydrogenated derivatives or 
    various substitutions.

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
    
    # Refined SMARTS patterns for beta-carboline core structure
    # Including possible hydrogenation and core variations
    beta_carboline_patterns = [
        'n1c2ccccc2c3[nH]c1c3',  # Fully aromatic core
        'Cn1c2ccccc2c3nccc3c1',  # Methylated indole nitrogen
        'n1c2ccccc2c3ccc[nH]c13',  # Another variation in hydrogenation
        'n1c2ccccc2c3[nH]ccc13',  # Original pattern
        '[nH]1c2ccccc2c3cccnc13',  # Another possible pattern with heterocycles
    ]

    # Check for any beta-carboline pattern match
    for pattern in beta_carboline_patterns:
        matcher = Chem.MolFromSmarts(pattern)
        if mol.HasSubstructMatch(matcher):
            return True, f"Contains a beta-carboline structure matching the pattern: {pattern}"

    # If no pattern matches
    return False, "No beta-carboline structure found matching the refined patterns"

# Example usage with one of the provided beta-carboline SMILES
smiles_example = "CCCNC(=O)N1CC2(C1)CN([C@H](C3=C2C4=C(N3C)C=C(C=C4)OC)CO)S(=O)(=O)C"
print(is_beta_carbolines(smiles_example))