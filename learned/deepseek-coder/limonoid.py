"""
Classifies: CHEBI:39434 limonoid
"""
"""
Classifies: CHEBI:26154 limonoid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_limonoid(smiles: str):
    """
    Determines if a molecule is a limonoid based on its SMILES string.
    Limonoids are highly oxygenated triterpenoids with a 4,4,8-trimethyl-17-furanylsteroid skeleton.

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

    # Check for furan ring (c1ccoc1)
    furan_pattern = Chem.MolFromSmarts("c1ccoc1")
    if not mol.HasSubstructMatch(furan_pattern):
        return False, "No furan ring found"

    # Check oxygen count (limonoids are highly oxygenated)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 5:
        return False, f"Too few oxygens ({o_count}), limonoids are highly oxygenated"

    # Check carbon count (triterpenoids have ~30 carbons)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 25 or c_count > 35:
        return False, f"Carbon count ({c_count}) not in triterpenoid range"

    # Check for 4,4,8-trimethyl pattern
    trimethyl_pattern = Chem.MolFromSmarts("[C;H3][C;H3][C;H3]")
    if not mol.HasSubstructMatch(trimethyl_pattern):
        return False, "No 4,4,8-trimethyl pattern found"

    # Check for steroid-like structure (4 fused rings)
    ring_info = mol.GetRingInfo()
    if len(ring_info.AtomRings()) < 4:
        return False, "Not enough rings for steroid-like structure"

    # Check molecular weight (should be >400 for limonoids)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f}) for limonoid"

    return True, "Contains furan ring, high oxygen content, triterpenoid structure, and 4,4,8-trimethyl pattern"