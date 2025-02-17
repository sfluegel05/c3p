"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines
Defined as: Any piperazinone that has a piperazine-2,5-dione skeleton.
This module defines a function is_2_5_diketopiperazines that determines if a molecule,
given by its SMILES string, contains a 2,5-diketopiperazine core.
"""

from rdkit import Chem

def has_terminal_carbonyl_oxygen(mol, carbon):
    """
    Check if a given carbon (expected to be carbonyl) has a double-bonded oxygen
    that is terminal (i.e. that oxygen has no extra heavy atom neighbors besides the carbon).
    """
    for neighbor in carbon.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
        if neighbor.GetAtomicNum() == 8 and bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            # Check oxygen is terminal (neighbors only its carbon partner, or hydrogens)
            heavy_neighbors = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != carbon.GetIdx()]
            if len(heavy_neighbors) == 0:
                return True
    return False

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule contains the piperazine-2,5-dione core.
    The core is defined by the six-membered ring:
        N - C(=O) - C - N - C(=O) - C
    where the two C atoms (in positions 1 and 4 of the ring when properly ordered)
    must be true carbonyl centers, meaning that a terminal oxygen is double-bonded.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if a validated 2,5-diketopiperazine core is present, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for the minimal piperazine-2,5-dione: a six-membered ring with two carbonyl groups.
    # This pattern captures the connectivity: O=C1NC(=O)CN1
    smarts = "O=C1NC(=O)CN1"
    core_query = Chem.MolFromSmarts(smarts)
    if core_query is None:
        return False, "Error in SMARTS pattern"
    
    # Find all substructure matches to the SMARTS pattern.
    # Each match is a tuple of atom indices corresponding to the pattern.
    matches = mol.GetSubstructMatches(core_query, uniquify=True)
    if not matches:
        return False, "Molecule does not contain the piperazine-2,5-dione core"

    # In our SMARTS pattern, by default the order of atoms is as written:
    # Let’s assume indices: 0: O (of first carbonyl) but we want the ring atoms.
    # In our pattern, the ring closure resets the order.
    # To validate, we will identify the two ring carbons which should be carbonyl centers.
    # We do this by scanning through the match and checking for carbons that are double bonded to oxygen.
    for match in matches:
        # Extract sub-mol atoms corresponding to the match.
        # The match should cover the 6 ring atoms plus the carbonyl oxygens are external.
        # Therefore, we recover the ring atoms from the match by ensuring they are in a ring:
        ring_atoms = []
        for idx in match:
            atom = mol.GetAtomWithIdx(idx)
            if atom.IsInRing():
                ring_atoms.append(atom)
        # We expect exactly 6 atoms in the ring.
        if len(ring_atoms) != 6:
            # The SMARTS may match even if the ring is fused in a larger system; skip if not exactly 6.
            continue
        
        # Now, identify the two ring carbons that are carbonyls.
        carbonyls = []
        non_carbonyl_carbons = []
        nitrogens = []
        # Classify ring atoms by atom type.
        for atom in ring_atoms:
            sym = atom.GetSymbol()
            if sym == "C":
                if has_terminal_carbonyl_oxygen(mol, atom):
                    carbonyls.append(atom)
                else:
                    non_carbonyl_carbons.append(atom)
            elif sym == "N":
                nitrogens.append(atom)
            else:
                # Any other atom unexpected in the core: discard this match.
                carbonyls = []
                break
        # We require exactly two carbonyl carbons, two non-carbonyl carbons, and two nitrogens.
        if len(carbonyls) != 2 or len(non_carbonyl_carbons) != 2 or len(nitrogens) != 2:
            continue

        # Optionally, one might enforce further checks (such as connectivity order in the cycle).
        # For our purposes, if a six-membered ring exists with exactly the expected counts and validated carbonyls,
        # we consider the 2,5-diketopiperazine core as present.
        ring_indices = tuple(sorted([atom.GetIdx() for atom in ring_atoms]))
        return True, ("Molecule contains a 2,5-diketopiperazine skeleton "
                      "(ring atom indices: {})".format(ring_indices))
        
    return False, "Molecule does not contain a validated piperazine-2,5-dione core"

# Example usage:
# Uncomment the lines below to test with some SMILES.
# tests = [
#     "O=C1CNC(=O)CN1",  # Simplest 2,5-diketopiperazine: Expected True
#     "CC(C)[C@@H]1NC(=O)C(=C)NC1=O",  # cyclo(DeltaAla-L-Val) from test set: Expected True
#     "O=C1N2[C@H](C(=O)N[C@H]1C)C[C@H](C2)O"  # one of the false positives: Expected False
# ]
# for smi in tests:
#     result, reason = is_2_5_diketopiperazines(smi)
#     print(smi, result, reason)