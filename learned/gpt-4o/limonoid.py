"""
Classifies: CHEBI:39434 limonoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids containing a 4,4,8-trimethyl-17-furanylsteroid skeleton.

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
    
    # Look for presence of a furan ring
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"
    
    # Check for high oxygenation; for simplification, let's assume at least 7 oxygens
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 7:
        return False, f"Insufficient oxygenation, found {o_count} oxygens"
    
    # Check if it's a triterpenoid by considering a potentially large number of carbon atoms (25-30)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (25 <= c_count <= 30):
        return False, f"Not a typical triterpenoid carbon count, found {c_count} carbons"
    
    # An additional check could be for a polycyclic framework
    n_rings = rdMolDescriptors.CalcNumRings(mol)
    if n_rings < 4:
        return False, f"Too few rings, found {n_rings}"
    
    return True, "Contains features typical of a limonoid (furan, high oxygenation, and polycyclic structure)"