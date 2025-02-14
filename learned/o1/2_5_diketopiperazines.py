"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines
"""

from rdkit import Chem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is defined as any piperazinone that has a piperazine-2,5-dione skeleton,
    which is a six-membered ring with nitrogen atoms at positions 1 and 4, and carbonyl groups at
    positions 2 and 5.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 2,5-diketopiperazine, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the 2,5-diketopiperazine core
    # Six-membered ring with two nitrogen atoms at positions 1 and 4,
    # carbonyl groups at positions 2 and 5 adjacent to the nitrogens
    diketopiperazine_pattern = Chem.MolFromSmarts('C1(=O)NC[C;R][N;R]C1=O')
    if diketopiperazine_pattern is None:
        return None, "Invalid SMARTS pattern for diketopiperazine"

    # Find all matching substructures in the molecule
    matches = mol.GetSubstructMatches(diketopiperazine_pattern)

    if matches:
        # Additional checks to confirm the presence of the correct ring system
        for match in matches:
            atoms_in_ring = [mol.GetAtomWithIdx(idx) for idx in match]

            # Get ring information
            ri = mol.GetRingInfo()
            ring_atoms = set()
            for ring in ri.AtomRings():
                if all(idx in ring for idx in match):
                    ring_atoms.update(ring)
                    break

            # Check that the ring is six-membered
            if len(ring_atoms) != 6:
                continue

            # Count nitrogen and carbon atoms in the ring
            n_count = 0
            c_count = 0
            for idx in ring_atoms:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 7:
                    n_count += 1
                elif atom.GetAtomicNum() == 6:
                    c_count += 1

            if n_count != 2 or c_count != 4:
                continue

            # Check that the carbonyl groups are at the correct positions
            carbonyl_positions = []
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6 and atom.IsInRing():
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8 and atom.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondTypeAsDouble() == 2.0:
                            carbonyl_positions.append(atom.GetIdx())

            if len(carbonyl_positions) != 2:
                continue

            # If all checks passed
            return True, "Molecule contains piperazine-2,5-dione skeleton"
        
        # If no matching rings after additional checks
        return False, "Molecule does not contain piperazine-2,5-dione skeleton"
    else:
        return False, "Molecule does not contain piperazine-2,5-dione skeleton"