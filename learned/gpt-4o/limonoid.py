"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids potentially containing a 4,4,8-trimethyl-17-furanylsteroid skeleton.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a limonoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for presence of a furan-like ring structure
    furan_patterns = [
        Chem.MolFromSmarts("c1ccoc1"),  # traditional furan
        Chem.MolFromSmarts("C1=COC=C1"),  # open to variations
    ]
    furan_found = any(mol.HasSubstructMatch(fp) for fp in furan_patterns)
    if not furan_found:
        return False, "No furan-like ring found"
    
    # Check for high oxygenation; previous rigid limit was inadequate
    # Check for at least 6 oxygens, but adapt to context
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 6:
        return False, f"Insufficient oxygenation, found {o_count} oxygens"
    
    # Evaluate carbon count more flexibly taking potential extensions
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (25 <= c_count <= 40):
        return False, f"Not a typical triterpenoid carbon count, found {c_count} carbons"
    
    # Examine adequate complexity for polycyclic framework
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, f"Too few rings, found {n_rings}"
    
    additional_pattern = Chem.MolFromSmarts("C1(C)C(C)(C)[C@H]2C[C@H](...)[C@@H]3")  # placeholder for triterpenoid backbone mimicry
    if not mol.HasSubstructMatch(additional_pattern):
        return False, "Missing typical backbone structure patterns"
    
    return True, "Contains features typical of a limonoid (furan-like structure, high oxygenation, and polycyclic complexity)"