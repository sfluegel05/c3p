"""
Classifies: CHEBI:72600 spiroketal
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_spiroketal(smiles: str):
    """
    Determines if a molecule is a spiroketal based on its SMILES string.
    A spiroketal is a cyclic ketal in which the ketal carbon is the only common atom of two rings.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a spiroketal, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find spiro atoms
    spiro_atoms = Chem.FindSpiroAtoms(mol)
    if not spiro_atoms:
        return False, "No spiro centers found"

    # Loop over spiro atoms to check for ketal functionality
    for atom_idx in spiro_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)

        # Check if the atom is carbon
        if atom.GetAtomicNum() != 6:
            continue  # Not a carbon atom

        # Get neighbors of the atom
        neighbors = atom.GetNeighbors()

        # Count oxygen atoms connected to this atom
        oxygens = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
        if len(oxygens) != 2:
            continue  # Not connected to exactly two oxygens

        # Check that the oxygens are connected to carbons (ether linkages)
        is_ketal = True
        for oxygen in oxygens:
            # Oxygen should be connected to two atoms: the spiro carbon and another carbon
            if oxygen.GetDegree() != 2:
                is_ketal = False
                break
            # Get the other atom connected to oxygen (excluding spiro carbon)
            oxygen_neighbors = [nbr for nbr in oxygen.GetNeighbors() if nbr.GetIdx() != atom_idx]
            if len(oxygen_neighbors) != 1:
                is_ketal = False
                break
            if oxygen_neighbors[0].GetAtomicNum() != 6:
                is_ketal = False
                break
        if not is_ketal:
            continue  # Does not fulfill ketal criteria

        # Check that spiro carbon has four single bonds and no hydrogens
        if atom.GetTotalDegree() != 4 or atom.GetTotalNumHs() != 0:
            continue  # Not a tetrahedral carbon without hydrogens

        # If all checks passed, this is a spiroketal
        return True, f"Molecule is a spiroketal with spiro ketal center at atom index {atom_idx}"

    return False, "No spiroketal functionality found"