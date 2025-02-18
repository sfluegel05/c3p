"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines (Any piperazinone that has a piperazine-2,5-dione skeleton.)

This module defines a function is_2_5_diketopiperazines that, given a SMILES string,
returns a boolean and an explanation string. It determines whether the molecule contains
a 2,5-diketopiperazine core defined as:
    N - C(=O) - C - N - C(=O) - C
arranged in a six-membered ring.
"""

from rdkit import Chem

def has_carbonyl(atom):
    """
    Check if the given carbon atom is part of a carbonyl group.
    That is, does it have at least one neighbor oxygen with a double bond?
    """
    if atom.GetSymbol() != "C":
        return False
    for neighbor in atom.GetNeighbors():
        bond = atom.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
        if neighbor.GetSymbol() == "O" and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            return True
    return False

def check_ring_pattern(atoms):
    """
    Given a list of 6 atoms (in order around the ring), check if they match the pattern:
    [N, C(=O), C(non-carbonyl), N, C(=O), C(non-carbonyl)]
    """
    # positions:
    # pos0: should be N
    if atoms[0].GetSymbol() != "N":
        return False
    # pos1: should be carbonyl carbon
    if not (atoms[1].GetSymbol() == "C" and has_carbonyl(atoms[1])):
        return False
    # pos2: should be carbon but not a carbonyl (it might be substituted though)
    if atoms[2].GetSymbol() != "C" or has_carbonyl(atoms[2]):
        return False
    # pos3: should be N
    if atoms[3].GetSymbol() != "N":
        return False
    # pos4: should be carbonyl carbon
    if not (atoms[4].GetSymbol() == "C" and has_carbonyl(atoms[4])):
        return False
    # pos5: should be carbon but not carbonyl
    if atoms[5].GetSymbol() != "C" or has_carbonyl(atoms[5]):
        return False
    return True

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is defined as any piperazinone that has a piperazine-2,5-dione
    skeleton. The core structure consists of a six-membered ring with the connectivity:
        N - C(=O) - C - N - C(=O) - C
    where the carbonyl groups (C with a =O) occur at positions 2 and 5 in the ring.
    
    This implementation checks all six-membered rings in the molecule and tests all cyclic
    orderings (rotations and reverse orderings) to see if one matches the above pattern.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule has a 2,5-diketopiperazine skeleton, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # returns tuples of atom indices
    
    # Look for a six-membered ring that matches our pattern.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # not six-membered, skip
        # Retrieve the corresponding atom objects.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Because the ring is cyclic, try all rotations and both directions:
        for start in range(6):
            # Extract ordering in forward direction:
            forward = [ring_atoms[(start + i) % 6] for i in range(6)]
            if check_ring_pattern(forward):
                return True, "Molecule contains a 2,5-diketopiperazine skeleton (found in ring with indices {})".format(ring)
            # Extract ordering in reverse direction:
            reverse = [ring_atoms[(start - i) % 6] for i in range(6)]
            if check_ring_pattern(reverse):
                return True, "Molecule contains a 2,5-diketopiperazine skeleton (found in ring with indices {})".format(ring)
    
    return False, "Molecule does not contain the required piperazine-2,5-dione skeleton"

# Example usage:
# Uncomment the following lines to test some examples:
# test_smiles = [
#     "O=C1CNC(=O)CN1",  # simplest 2,5-diketopiperazine, should be True
#     "CC1C(=O)NC(=O)N1",  # 5-membered ring – should be False
# ]
# for smi in test_smiles:
#     result, reason = is_2_5_diketopiperazines(smi)
#     print(smi, result, reason)