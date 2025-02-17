"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines (any piperazinone that has a piperazine-2,5-dione skeleton)

This module defines a function is_2_5_diketopiperazines that determines whether a molecule,
given by its SMILES string, contains a 2,5-diketopiperazine core. The core is defined by the pattern:
    N - C(=O) - C - N - C(=O) - C
arranged as a six-membered ring.
"""

from rdkit import Chem

def has_carbonyl(mol, atom):
    """
    Check if the specified carbon atom is part of a carbonyl group.
    That is, does it have at least one neighboring oxygen connected by a double bond?

    Args:
        mol (rdkit.Chem.Mol): The molecule to which the atom belongs.
        atom (rdkit.Chem.Atom): The atom to check.

    Returns:
        bool: True if a double-bonded oxygen is found, False otherwise.
    """
    if atom.GetSymbol() != "C":
        return False
    # Loop through all neighbors of the atom
    for neighbor in atom.GetNeighbors():
        # Get the bond between the two atoms via the molecule
        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
        if neighbor.GetSymbol() == "O" and bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            return True
    return False

def check_ring_pattern(mol, atoms):
    """
    Given a list of 6 atoms in order around a ring, check if they match the pattern:
    [N, C(=O), C (non-carbonyl), N, C(=O), C (non-carbonyl)].

    Args:
        mol (rdkit.Chem.Mol): The molecule.
        atoms (list of rdkit.Chem.Atom): List of atoms in the ring (ordered).

    Returns:
        bool: True if the atoms match the desired pattern, False otherwise.
    """
    # Position 0: Must be nitrogen
    if atoms[0].GetSymbol() != "N":
        return False
    # Position 1: Must be a carbonyl carbon (C with a double-bonded oxygen)
    if not (atoms[1].GetSymbol() == "C" and has_carbonyl(mol, atoms[1])):
        return False
    # Position 2: Must be a carbon that is not a carbonyl
    if atoms[2].GetSymbol() != "C" or has_carbonyl(mol, atoms[2]):
        return False
    # Position 3: Must be nitrogen
    if atoms[3].GetSymbol() != "N":
        return False
    # Position 4: Must be a carbonyl carbon
    if not (atoms[4].GetSymbol() == "C" and has_carbonyl(mol, atoms[4])):
        return False
    # Position 5: Must be a carbon that is not a carbonyl
    if atoms[5].GetSymbol() != "C" or has_carbonyl(mol, atoms[5]):
        return False
    return True

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine has a six-membered ring with the connectivity:
        N - C(=O) - C - N - C(=O) - C,
    where the two carbons in positions 1 and 4 are carbonyl centers (i.e., have a double-bonded oxygen).

    The function checks all six-membered rings in the molecule and all rotations (and their reversals)
    to see if one of them matches the specified pattern.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule has a 2,5-diketopiperazine skeleton, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()  # Tuples of atom indices forming rings.
    
    # Check each six-membered ring for the 2,5-diketopiperazine pattern.
    for ring in atom_rings:
        if len(ring) != 6:
            continue  # Only interested in six-membered rings.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Because the ring is cyclic, test all cyclic rotations and their reverse orders.
        for start in range(6):
            # Forward ordering
            forward = [ring_atoms[(start + i) % 6] for i in range(6)]
            if check_ring_pattern(mol, forward):
                return True, "Molecule contains a 2,5-diketopiperazine skeleton (ring indices: {})".format(ring)
            # Reverse ordering (to account for opposite direction around the ring)
            reverse = [ring_atoms[(start - i) % 6] for i in range(6)]
            if check_ring_pattern(mol, reverse):
                return True, "Molecule contains a 2,5-diketopiperazine skeleton (ring indices: {})".format(ring)
    
    return False, "Molecule does not contain the required piperazine-2,5-dione skeleton"

# Example usage (you can uncomment these lines to test a couple of cases):
# test_smiles = [
#     "O=C1CNC(=O)CN1",  # Simplest 2,5-diketopiperazine, expected: True
#     "CC1C(=O)NC(=O)N1",  # Not a six-membered ring as required, expected: False
# ]
# for smi in test_smiles:
#     result, reason = is_2_5_diketopiperazines(smi)
#     print(smi, result, reason)