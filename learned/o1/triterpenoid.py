"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:26870 triterpenoid
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is derived from a triterpene, typically containing 30 carbons,
    but may have modifications such as rearrangements, glycosylations, or additions/removals of methyl groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpenoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Exclude molecules that are steroids
    steroid_smarts = "[#6]12[#6][#6][#6]3[#6]([#6]1)[#6]4[#6][#6][#6][#6][#6]4[#6]3[#6][#6]2"
    steroid_pattern = Chem.MolFromSmarts(steroid_smarts)
    if mol.HasSubstructMatch(steroid_pattern):
        return False, "Molecule contains steroid nucleus, possibly a steroid, not a triterpenoid"
    
    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20 or c_count > 60:
        return False, f"Carbon count {c_count} not in typical triterpenoid range (20-60 carbons)"
    
    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, less likely to be a triterpenoid"
    
    # Check for presence of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings == 0:
        return False, "No rings found, less likely to be a triterpenoid"
    
    # Optionally, check for isoprene units
    # Triterpenoids are formed from six isoprene units (C5H8)
    # This is an advanced check and may not be reliable for all structures
    
    # If the molecule passes the above checks, classify as triterpenoid
    return True, "Molecule meets criteria for a triterpenoid (carbon count, oxygen atoms, ring structures)"