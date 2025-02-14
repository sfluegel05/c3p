"""
Classifies: CHEBI:65061 2,5-diketopiperazines
"""
"""
Classifies: 2,5-diketopiperazines
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_2_5_diketopiperazines(smiles: str):
    """
    Determines if a molecule is a 2,5-diketopiperazine based on its SMILES string.
    A 2,5-diketopiperazine is any piperazinone that has a piperazine-2,5-dione skeleton,
    which is a six-membered ring with nitrogen atoms at positions 1 and 4, and carbonyl
    groups at positions 2 and 5.

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

    # Get ring information
    ri = mol.GetRingInfo()
    rings = ri.AtomRings()

    # Iterate over all rings in the molecule
    for ring in rings:
        if len(ring) != 6:
            continue  # Skip if not a six-membered ring

        # Initialize counters for nitrogen and carbon atoms
        n_nitrogens = 0
        n_carbons = 0
        nitrogen_atoms = []
        carbon_atoms = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 7:
                n_nitrogens += 1
                nitrogen_atoms.append(idx)
            elif atomic_num == 6:
                n_carbons += 1
                carbon_atoms.append(idx)
            else:
                # Ring contains atoms other than carbon or nitrogen
                break
        else:
            # Only proceed if the ring has exactly 2 nitrogens and 4 carbons
            if n_nitrogens != 2 or n_carbons != 4:
                continue

            # Check that each nitrogen atom is adjacent to a carbonyl carbon
            has_carbonyls = True
            for n_idx in nitrogen_atoms:
                n_atom = mol.GetAtomWithIdx(n_idx)
                neighbor_in_ring = False
                for neighbor in n_atom.GetNeighbors():
                    if neighbor.GetIdx() in ring and neighbor.GetAtomicNum() == 6:
                        neighbor_in_ring = True
                        c_atom = neighbor
                        # Check if carbon has a double bond to oxygen
                        has_carbonyl = False
                        for c_neighbor in c_atom.GetNeighbors():
                            bond = mol.GetBondBetweenAtoms(c_atom.GetIdx(), c_neighbor.GetIdx())
                            if c_neighbor.GetAtomicNum() == 8 and bond.GetBondType() == rdchem.BondType.DOUBLE:
                                has_carbonyl = True
                                break
                        if not has_carbonyl:
                            has_carbonyls = False
                            break
                if not neighbor_in_ring or not has_carbonyls:
                    has_carbonyls = False
                    break

            if has_carbonyls:
                return True, "Molecule contains piperazine-2,5-dione skeleton"

    return False, "Molecule does not contain piperazine-2,5-dione skeleton"