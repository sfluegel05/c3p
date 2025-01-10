"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are terpenoids consisting of 25 carbon atoms (five isoprene units),
    often with complex ring structures and various functional groups.

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
    
    # Count carbon atoms to check if it has around 25 carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Allow a flexible range to cover edge cases, e.g., 18-32
    if not (18 <= c_count <= 32):
        return False, "Does not have a typical carbon count for sesterterpenoids"
    
    # Ensure the molecule contains multiple ring structures
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings < 3:
        return False, "Not enough ring structures; typically sesterterpenoids have multiple rings"
    
    # Hypothetical SMARTS pattern to find common terpenoid-like structures
    # This is a placeholder and should be replaced with actual SMARTS based on known structures
    terpenoid_pattern = Chem.MolFromSmarts("C1CCCC1") # Simplistic example
    if not mol.HasSubstructMatch(terpenoid_pattern):
        return False, "No common terpenoid structure found"
   
    # Check for typical features (e.g., oxygen presence in rings, stereochemistry)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found, typically present in sesterterpenoids"
    
    # If carbon count, ring structures and additional checks pass, posit it as a sesterterpenoid
    return True, "Carbon count, ring structures, and features suggest a sesterterpenoid"