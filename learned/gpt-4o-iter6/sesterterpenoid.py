"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 25:
        return False, f"Carbon count ({c_count}) not typical for a sesterterpenoid (usually around 25)"

    # Check for typical terpenoid structure features
    # Terpenoids are often composed of isoprene units (C5H8), seeking patterns
    isoprene_pattern = Chem.MolFromSmarts("C=C(C)C")
    if not mol.HasSubstructMatch(isoprene_pattern):
        return False, "No isoprene units found which are common in terpenoids"

    # Additional checks can be implemented as necessary
    # For complex cases, consider manual curation and expert judgment
    
    # If C25 skeleton and terpenoid features present, classify as sesterterpenoid
    return True, "Likely a sesterterpenoid based on carbon count and structural features"