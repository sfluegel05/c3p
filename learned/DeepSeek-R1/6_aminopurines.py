"""
Classifies: CHEBI:20706 6-aminopurines
"""
"""
Classifies: CHEBI:16761 6-aminopurine derivatives
"""
from rdkit import Chem

def is_6_aminopurines(smiles: str):
    """
    Determines if a molecule is a 6-aminopurine derivative based on its SMILES string.
    A 6-aminopurine contains the adenine substructure (6-aminopurine core).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains adenine substructure, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define adenine (6-aminopurine) substructure pattern with explicit 6-amino group
    # SMARTS pattern for adenine: [NH2] attached to position 6 of purine ring
    adenine_pattern = Chem.MolFromSmarts("[NH2]c1nc2c(ncnc2n1)")
    if adenine_pattern is None:
        return False, "Invalid adenine pattern"

    # Check for presence of adenine substructure
    if mol.HasSubstructMatch(adenine_pattern):
        return True, "Contains 6-aminopurine (adenine) substructure"
    else:
        return False, "No adenine substructure found"