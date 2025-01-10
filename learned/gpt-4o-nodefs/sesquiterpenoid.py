"""
Classifies: CHEBI:26658 sesquiterpenoid
"""
from rdkit import Chem

def is_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a sesquiterpenoid based on its SMILES string.
    A sesquiterpenoid typically has a 15-carbon skeleton derived from three isoprene units,
    often arranged in complex ring structures.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesquiterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for approximately 15 carbons, allowing some variance
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (13 <= c_count <= 17):
        return False, f"Carbon count is {c_count}, expected closer to 15 for sesquiterpenoid"

    # Identify sesquiterpenoid-specific ring patterns
    # Example pattern: bicyclic structures (proxy)
    bicyclic_pattern = Chem.MolFromSmarts("C1CCC2CCC(C1)C2")  # Example tweak
    if not mol.HasSubstructMatch(bicyclic_pattern):
        return False, "No complex bicyclic-like structure found"

    # Check for functional groups commonly found in terpenoids
    alcohol_pattern = Chem.MolFromSmarts("[CX4][OH]")
    lactone_pattern = Chem.MolFromSmarts("O=C1OC(C=CC1)C")
    
    if not mol.HasSubstructMatch(alcohol_pattern) and not mol.HasSubstructMatch(lactone_pattern):
        return False, "Missing common functional groups like alcohols or lactones"
    
    return True, "Contains 15-carbon skeleton with sesquiterpenoid-like complex ring structures and functional groups"