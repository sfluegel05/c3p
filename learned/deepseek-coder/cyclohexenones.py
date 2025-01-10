"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: CHEBI:36423 cyclohexenone
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is a six-membered alicyclic ketone with one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule has a six-membered ring
    ring_info = mol.GetRingInfo()
    six_membered_rings = [ring for ring in ring_info.AtomRings() if len(ring) == 6]
    if not six_membered_rings:
        return False, "No six-membered ring found"

    # Check for the presence of a ketone group (C=O) in the ring
    ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone group found"

    # Check if the ketone is part of the six-membered ring
    ketone_in_ring = any(any(atom_idx in ring for atom_idx in match) for match in ketone_matches for ring in six_membered_rings)
    if not ketone_in_ring:
        return False, "Ketone group not in the six-membered ring"

    # Check for exactly one double bond in the ring
    double_bond_count = 0
    for ring in six_membered_rings:
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring and bond.GetBondType() == Chem.BondType.DOUBLE:
                double_bond_count += 1
    if double_bond_count != 1:
        return False, f"Found {double_bond_count} double bonds in the ring, need exactly 1"

    return True, "Contains a six-membered ring with one double bond and a ketone group"