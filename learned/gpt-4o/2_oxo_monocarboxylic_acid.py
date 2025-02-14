"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.

    A 2-oxo monocarboxylic acid has:
    - A carbonyl group (C=O) at the alpha position.
    - A carboxylic acid functional group (C(=O)O).
    - Any R-group can be attached to the carbon bearing the oxo group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS pattern for 2-oxo monocarboxylic acids
    patterns = [
        Chem.MolFromSmarts("C(=O)[C][C](=O)O"),  # General form
        Chem.MolFromSmarts("C(=O)[CH2][C](=O)O"),  # CH2 next to C=O
        Chem.MolFromSmarts("[#6](=O)[C][C](=O)O")  # More generalized form with R
    ]

    # Check if any of the patterns match
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a 2-oxo monocarboxylic acid substructure"

    return False, "Does not contain a 2-oxo monocarboxylic acid substructure"

# Example usage
# smiles = "CCC(=O)C(O)=O"  # 2-oxobutanoic acid
# print(is_2_oxo_monocarboxylic_acid(smiles))