"""
Classifies: CHEBI:13601 3-oxo-5alpha-steroid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_3_oxo_5alpha_steroid(smiles: str):
    """
    Determines if a molecule is a 3-oxo-5alpha-steroid based on its SMILES string.
    A 3-oxo-5alpha-steroid has a steroid nucleus, a ketone at position 3, and alpha configuration at position 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-oxo-5alpha-steroid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of a ketone group (C=O)
    ketone_pattern = Chem.MolFromSmarts("[CX3]=[OX1]")
    if not mol.HasSubstructMatch(ketone_pattern):
        return False, "No ketone group found"

    # Check for steroid nucleus (three 6-membered rings and one 5-membered ring)
    ring_info = mol.GetRingInfo()
    ring_sizes = [len(ring) for ring in ring_info.AtomRings()]
    sorted_ring_sizes = sorted(ring_sizes)
    if sorted_ring_sizes.count(6) < 3 or 5 not in sorted_ring_sizes:
        return False, "Incorrect ring structure for steroid nucleus"

    # Check for 5alpha configuration (trans-decalin pattern in A/B rings)
    # SMARTS pattern for trans-decalin (5alpha configuration)
    trans_decalin_pattern = Chem.MolFromSmarts("[C@]1([C@@]2(CCCC[C@H]2C)C)[C@H]1C")
    if not mol.HasSubstructMatch(trans_decalin_pattern):
        return False, "A/B ring fusion not in 5alpha configuration"

    # Verify ketone is in a six-membered ring (A ring)
    ketone_in_A_ring = False
    for atom_idx in mol.GetSubstructMatches(ketone_pattern):
        atom = mol.GetAtomWithIdx(atom_idx[0])
        for ring in ring_info.AtomRings():
            if len(ring) == 6 and atom_idx[0] in ring:
                ketone_in_A_ring = True
                break
        if ketone_in_A_ring:
            break
    if not ketone_in_A_ring:
        return False, "Ketone not located in a six-membered ring (A ring)"

    return True, "3-oxo group and 5alpha steroid configuration confirmed"