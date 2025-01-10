"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is a compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for primary amine: -NH2 attached to a carbon
    primary_amine_pattern = Chem.MolFromSmarts("[NH2][CX4]")
    
    # Check if the molecule matches the primary amine pattern
    if mol.HasSubstructMatch(primary_amine_pattern):
        return True, "Contains a primary amine group (-NH2 attached to a carbon)"
    else:
        return False, "No primary amine group (-NH2 attached to a carbon) found"