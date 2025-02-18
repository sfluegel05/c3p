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

def has_terminal_carbonyl(mol, carbon):
    """
    Checks if the given carbon atom is part of a carbonyl that has a terminal oxygen.
    Terminal means that the oxygen atom double-bonded to the carbon has no additional heavy atom neighbors.
    """
    # Loop over neighbors of the carbon
    for neighbor in carbon.GetNeighbors():
        bond = mol.GetBondBetweenAtoms(carbon.GetIdx(), neighbor.GetIdx())
        # Check for a double-bonded oxygen atom
        if neighbor.GetAtomicNum() == 8 and bond and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            # Verify that the oxygen's other connections are only to hydrogens (or none)
            heavy_neighbors = [
                nbr for nbr in neighbor.GetNeighbors() 
                if nbr.GetAtomicNum() > 1 and nbr.GetIdx() != carbon.GetIdx()
            ]
            if len(heavy_neighbors) == 0:
                return True
    return False

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule contains a 2,5-diketopiperazine core based on its SMILES.
    The 2,5-diketopiperazine core is defined as a six-membered ring with a pattern
    (read in cyclic order): N, C(carbonyl), C(non-carbonyl), N, C(carbonyl), C(non-carbonyl),
    where the two carbonyl carbons have a terminal oxygen double-bond.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if a 2,5-diketopiperazine core is detected, False otherwise.
        str: Explanation for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings using RDKit's symmetric SSSR algorithm.
    rings = Chem.GetSymmSSSR(mol)
    
    # Process each ring of size 6.
    for ring in rings:
        if len(ring) != 6:
            continue  # Only interested in 6-membered rings
        # Convert ring (which is a tuple of atom indices) to a list of atoms.
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        
        # The desired pattern (in cyclic order) is:
        # [0] N, [1] C with terminal carbonyl, [2] C without carbonyl,
        # [3] N, [4] C with terminal carbonyl, [5] C without carbonyl.
        # As the ring can be rotated arbitrarily, we try all rotations.
        for shift in range(6):
            ordered = ring_atoms[shift:] + ring_atoms[:shift]
            valid = True
            # Position 0 must be a nitrogen.
            if ordered[0].GetSymbol() != "N":
                valid = False
            # Position 1 must be a carbon with a terminal C=O.
            if not (ordered[1].GetSymbol() == "C" and has_terminal_carbonyl(mol, ordered[1])):
                valid = False
            # Position 2 must be a carbon and should NOT itself be a carbonyl carbon.
            if ordered[2].GetSymbol() != "C" or has_terminal_carbonyl(mol, ordered[2]):
                valid = False
            # Position 3 must be a nitrogen.
            if ordered[3].GetSymbol() != "N":
                valid = False
            # Position 4 must be a carbon with a terminal C=O.
            if not (ordered[4].GetSymbol() == "C" and has_terminal_carbonyl(mol, ordered[4])):
                valid = False
            # Position 5 must be a carbon and not a carbonyl carbon.
            if ordered[5].GetSymbol() != "C" or has_terminal_carbonyl(mol, ordered[5]):
                valid = False

            if valid:
                ring_indices = [atom.GetIdx() for atom in ordered]
                return True, ("Molecule contains a 2,5-diketopiperazine core "
                              f"with ring atom indices (in order): {ring_indices}")
    
    return False, "Molecule does not contain a validated piperazine-2,5-dione core"

# Example usage:
# Uncomment the lines below to test with some SMILES.
# tests = [
#     "O=C1CNC(=O)CN1",  # simplest 2,5-diketopiperazine: Expected True
#     "CC(C)[C@@H]1NC(=O)C(=C)NC1=O",  # cyclo(DeltaAla-L-Val): Expected True (if the core is preserved)
#     "S1S[C@]23N([C@@H]4[C@@H](O)C=C[C@@H]([C@H]4C2)O)C([C@]15N([C@@H]6[C@@H](O)C=C[C@@H]([C@H]6C5)O)C3=O)=O"  # Brocazine E: Expected True if the core is present
# ]
# for smi in tests:
#     result, reason = is_2_5_diketopiperazines(smi)
#     print(smi, result, reason)