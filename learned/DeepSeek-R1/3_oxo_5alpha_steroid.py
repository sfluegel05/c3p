"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid has a steroid nucleus (three 6-membered and one 5-membered rings),
    a ketone at position 3, and alpha configuration (trans A/B ring fusion) at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of a ketone group (C=O) at position 3 in steroid skeleton
    # Position 3 is in the A ring (six-membered), adjacent to two carbons in the ring
    # SMARTS pattern for 3-oxo in steroid A ring: ketone connected to two adjacent carbons in a six-membered ring
    ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1;R]")
    ketone_matches = mol.GetSubstructMatches(ketone_pattern)
    if not ketone_matches:
        return False, "No ketone group found in ring system"

    # Verify ketone is in a six-membered ring (A ring)
    ring_info = mol.GetRingInfo()
    valid_ketone = False
    for atom_idx in ketone_matches:
        atom = mol.GetAtomWithIdx(atom_idx[0])
        for ring in ring_info.AtomRings():
            if len(ring) == 6 and atom_idx[0] in ring:
                # Check if ketone is at position 3 (adjacent to two adjacent carbons in the ring)
                neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
                if len(neighbors) >= 2:
                    # Check if neighbors are in the same ring and adjacent to each other
                    if (neighbors[0] in ring) and (neighbors[1] in ring):
                        valid_ketone = True
                        break
        if valid_ketone:
            break
    if not valid_ketone:
        return False, "Ketone not in correct position (C3 of A ring)"

    # Check steroid nucleus structure (three 6-membered rings and one 5-membered ring)
    ring_sizes = sorted([len(r) for r in ring_info.AtomRings()])
    if not (ring_sizes.count(6) >= 3 and 5 in ring_sizes):
        return False, "Incorrect ring structure for steroid nucleus"

    # Check 5alpha configuration (trans A/B ring fusion)
    # SMARTS pattern for trans-decalin system (5alpha configuration)
    trans_decalin_pattern = Chem.MolFromSmarts("[C@H]1[C@@H]2CCCC[C@H]2CC1")
    if not trans_decalin_pattern:
        return None, None  # Pattern creation failed
    if not mol.HasSubstructMatch(trans_decalin_pattern):
        return False, "A/B ring fusion not in 5alpha configuration"

    return True, "3-oxo group in A ring and 5alpha configuration confirmed"