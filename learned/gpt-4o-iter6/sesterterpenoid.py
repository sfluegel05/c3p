"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are derived from sesterterpenes and can have a modified C25 skeleton.

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

    # Count carbon atoms, consider modified C25 skeletons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 43:
        return False, f"Carbon count ({c_count}) not typical or too large for a sesterterpenoid"

    # Define broader patterns for isoprenoid-like structures, considering modifications
    # Typical sesterterpenoid may include complex ring systems
    isoprene_patterns = [
        Chem.MolFromSmarts("C=C(C)CC"),  # Basic isoprene unit
        Chem.MolFromSmarts("C=CC(C)C"),  # Variation of isoprene
        Chem.MolFromSmarts("CC(C)=C"),   # Further variation
        # Additional cyclic patterns could be checked here
    ]
    
    # Check if important cyclic structures indicative of sesterterpenoids are present
    cyclic_pattern = Chem.MolFromSmarts("C1CCC(CC1)C")  # Example cyclic pattern
    if not any(mol.HasSubstructMatch(pattern) for pattern in isoprene_patterns) and not mol.HasSubstructMatch(cyclic_pattern):
        return False, "No isoprene-like or key cyclic structures typical of sesterterpenoids found"

    # If structure fits broad criteria
    return True, "Likely a sesterterpenoid based on adjusted carbon count and structural features"