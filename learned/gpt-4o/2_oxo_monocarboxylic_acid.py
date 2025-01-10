"""
Classifies: CHEBI:35910 2-oxo monocarboxylic acid
"""
from rdkit import Chem

def is_2_oxo_monocarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-oxo monocarboxylic acid based on its SMILES string.
    A 2-oxo monocarboxylic acid has a 2-oxo substituent relative to a monocarboxylic acid group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2-oxo monocarboxylic acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Separate SMARTS patterns for oxo group and carboxylic acid group
    # Oxo group; consider both C(=O) attachments
    oxo_pattern = Chem.MolFromSmarts("C(=O)")
    # Monocarboxylic acid group C(=O)O
    carboxylic_pattern = Chem.MolFromSmarts("C(=O)O")

    # Match for both oxo and carboxylic acid group
    oxo_matches = mol.GetSubstructMatches(oxo_pattern)
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_pattern)
    
    if not oxo_matches or not carboxylic_matches:
        return False, "Missing essential group"

    # Check proximity and correct structure
    for oxo in oxo_matches:
        for carboxy in carboxylic_matches:
            # Ensure proximity in a flexible sense, oxo must be at C2 relative to carboxylic
            if abs(oxo[0] - carboxy[0]) in [1, 2]:  # Allow for moderate flexibility 
                return True, "Contains both 2-oxo group and monocarboxylic acid group"

    return False, "No 2-oxo substituent at correct location relative to monocarboxylic acid found"