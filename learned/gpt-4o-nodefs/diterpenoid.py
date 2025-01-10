"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid typically contains 20 carbon atoms and may exhibit a complex ring structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string using RDKit
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms
    num_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if num_carbons != 20:
        return False, f"Contains {num_carbons} carbon atoms, but diterpenoids should have exactly 20"
    
    # Ring structures are common in diterpenoids
    ring_info = mol.GetRingInfo()
    # Confirm presence of at least one large ring
    num_rings = ring_info.NumRings()
    if num_rings < 1:
        return False, "No ring structures found, diterpenoids typically have complex ring structures"
    
    # Check for the presence of functional groups (e.g., alcohol, carbonyl, epoxide)
    # Example: carbonyl group pattern
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    carbonyl_matches = mol.HasSubstructMatch(carbonyl_pattern)
    
    if not carbonyl_matches:
        return False, "Missing typical functional groups like carbonyl"
    
    return True, "Likely a diterpenoid based on carbon count and structural complexity"