"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids generally have a base skeleton of about 30 carbon atoms, possibly rearranged or modified.
    They often include multiple ring structures and functional groups like -OH, =O, or -COOH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms and validate the range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (25 <= c_count <= 35):
        return False, f"Carbon count {c_count} not in triterpenoid range (25-35)"
    
    # Check for presence of oxygen to imply functional groups
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; typical triterpenoids have oxygenated functional groups"
    
    # Ensure adequate functional groups are present
    if not any(atom.GetDegree() > 1 and atom.GetSymbol() == 'O' for atom in mol.GetAtoms()):
        return False, "No eligible functional groups indicative of triterpenoid modifications"

    # Verify presence of ring structures; triterpenoids usually have cyclic backbones
    ring_info = mol.GetRingInfo()
    if not ring_info.NumRings():
        return False, "No ring structures found; triterpenoids usually have rings"
    
    # If all checks pass, classify as a triterpenoid
    return True, "Molecule matches triterpenoid characteristics based on carbon count, functional groups, and ring presence"