"""
Classifies: CHEBI:36615 triterpenoid
"""
"""
Classifies: CHEBI:26870 triterpenoid
"""

from rdkit import Chem

def is_triterpenoid(smiles: str):
    """
    Determines if a molecule is a triterpenoid based on its SMILES string.
    A triterpenoid is derived from a triterpene, typically containing 30 carbons,
    but may have modifications such as rearrangements, glycosylations, or removal of methyl groups.

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

    # Exclude molecules containing atoms other than C, H, O
    allowed_atoms = {1, 6, 8}  # H, C, O
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Contains atoms other than C, H, O"

    # Count carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 18 or c_count > 70:
        return False, f"Carbon count {c_count} not in typical triterpenoid range (18-70 carbons)"

    # Count oxygen atoms
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count == 0:
        return False, "No oxygen atoms found, less likely to be a triterpenoid"
    if o_count > 20:
        return False, f"Too many oxygen atoms ({o_count}), may be a polysaccharide or not a triterpenoid"

    # Check for number of rings
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Only {num_rings} rings found, less likely to be a triterpenoid"

    # Check for aromatic rings (triterpenoids are usually non-aromatic)
    num_aromatic = sum(1 for ring in ring_info.AtomRings() if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring))
    if num_aromatic > 0:
        return False, "Contains aromatic rings, less likely to be a triterpenoid"

    # If the molecule passes the above checks, classify as triterpenoid
    return True, "Molecule meets criteria for a triterpenoid (atom types, carbon count, oxygen atoms, ring structures)"