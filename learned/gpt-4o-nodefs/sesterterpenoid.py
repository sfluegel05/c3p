"""
Classifies: CHEBI:26660 sesterterpenoid
"""
from rdkit import Chem

def is_sesterterpenoid(smiles: str):
    """
    Determines if a molecule is a sesterterpenoid based on its SMILES string.
    Sesterterpenoids are characterized by 25 carbon atoms, complex ring structures,
    and functional groups like alcohols, ketones, and olefins.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sesterterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
    
    # Carbon count typical for sesterterpenoids is around 25
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (24 <= c_count <= 26):
        return False, f"Carbon count {c_count} not typical for sesterterpenoids"
    
    # Ensure the molecule contains at least 3 rings
    num_rings = mol.GetRingInfo().NumRings()
    if num_rings < 3:
        return False, f"Not enough ring structures; found {num_rings} rings"
    
    # Check for sesterterpenoid-specific patterns
    sesterterpenoid_smarts = [
        Chem.MolFromSmarts("C1CCC2C(C1)CC(C)(C)C2"), # example of a common motif
        Chem.MolFromSmarts("C=CC=C(C)C=C"), # double bonds in carbon chains
        Chem.MolFromSmarts("[C&R](O)[C&R]"), # alcohol groups on rings
    ]
    
    if not any(mol.HasSubstructMatch(pattern) for pattern in sesterterpenoid_smarts):
        return False, "No characteristic sesterterpenoid patterns found"
    
    # Check for presence of heteroatoms like oxygen, which are typical in sesterterpenoids
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Insufficient presence of oxygen atoms, common in sesterterpenoids"

    # Given these checks, the molecule might be a sesterterpenoid
    return True, "Carbon, ring structures, and features suggest it's a sesterterpenoid"