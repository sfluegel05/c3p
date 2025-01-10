"""
Classifies: CHEBI:23849 diterpenoid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_diterpenoid(smiles: str):
    """
    Determines if a molecule is a diterpenoid based on its SMILES string.
    A diterpenoid has a C20 skeleton derived from a diterpene, potentially rearranged or modified.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a diterpenoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count the number of carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Check for approximate C20 skeleton
    if not (18 <= c_count <= 22):  # Allow for minor modifications
        return False, f"Carbon count ({c_count}) not typical for diterpenoid"

    # Check for terpenoid pattern (e.g., C5H8 isoprene units)
    isoprene_pattern = Chem.MolFromSmarts("[C]1([CH2])[CH2][CH][CH2]1")
    has_isoprene_units = mol.HasSubstructMatch(isoprene_pattern)

    if not has_isoprene_units:
        return False, "No isoprene-like units typical in terpenoids found"
    
    # Other pattern checks can be added here to further confirm

    return True, "Molecule matches diterpenoid characteristics with C20-like skeleton derived from diterpene"