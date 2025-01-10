"""
Classifies: CHEBI:36615 triterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    Triterpenoids generally have a base skeleton involving multiple rings and 
    distinctive functional groups such as hydroxyl, ketone, or carboxyl groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Count carbon atoms and allow a more lenient range
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 40:
        return False, f"Carbon count {c_count} is outside an extended range for typical triterpenoids (20-40)"
    
    # Check for a minimum number of oxygen atoms suggesting functional modifications
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 1:
        return False, "No oxygen atoms found; typical triterpenoids have oxygenated functional groups"

    # Identify specific functional groups often seen in triterpenoids
    if not any(atom.GetSymbol() == 'O' and atom.GetDegree() > 1 for atom in mol.GetAtoms()):
        return False, "No eligible functional groups indicative of triterpenoid modifications (e.g., -OH, =O)"

    # Verify the presence of significant ring structures
    ring_info = mol.GetRingInfo()
    if ring_info.NumRings() < 3:  # Triterpenoids often have multiple rings
        return False, f"Insufficient ring structures found, expected at least 3, found {ring_info.NumRings()}"

    # If all criteria are met, classify as a triterpenoid
    return True, "Molecule matches triterpenoid characteristics based on flexible carbon count, functional groups, and ring presence"