"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are terpenoids with about 25 carbon atoms, complex ring structures,
    and various functional groups.

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
    # Use a more precise carbon count range typical for sesterterpenoids
    if not (22 <= c_count <= 28):
        return False, f"Does not have typical carbon count for sesterterpenoids, found {c_count} carbons"
    
    # Ensure the molecule contains multiple ring structures
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings < 3:
        return False, f"Not enough ring structures; found {num_rings} rings"
    
    # More sophisticated terpenoid pattern including stereochemistry
    # Use a hypothetical SMARTS pattern more specific to terpenoids
    # Sesterterpenoids often have complex, interconnected rings
    terpenoid_patterns = [
        Chem.MolFromSmarts("[C&R]-[C&R]-[C&R](=O)-[C&R]-[C&R]"),  # Example pattern involving ketone and carbon atoms
        Chem.MolFromSmarts("C1CCC(CC1)C")  # Possible complex ring pattern
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in terpenoid_patterns):
        return False, "No terpenoid-like ring structures found"
    
    # Check for oxygen or other heteroatoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found, typically present in sesterterpenoids"
    
    # Given these checks, the molecule is hypothesized to be a sesterterpenoid
    return True, "Carbon count, ring structures, and features suggest a sesterterpenoid"