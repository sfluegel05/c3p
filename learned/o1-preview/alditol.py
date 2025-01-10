"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:18385 alditol
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH,
    formally derivable from an aldose by reduction of the carbonyl group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if the molecule is acyclic
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic"

    # Get all carbon atoms in the molecule
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) == 0:
        return False, "No carbon atoms found"

    # Find the main carbon chain (ignore branches for now)
    paths = Chem.rdmolops.FindAllPathsOfLengthN(mol, len(carbons), useBonds=False)
    if not paths:
        return False, "No continuous carbon chain found"

    # Choose the longest carbon chain
    max_path = max(paths, key=lambda p: len(p))
    chain_atoms = [mol.GetAtomWithIdx(idx) for idx in max_path]

    # Check terminal carbons for CH2OH groups
    termini = [chain_atoms[0], chain_atoms[-1]]
    for term_atom in termini:
        # Terminal carbon should be connected to CH2OH
        neighbors = term_atom.GetNeighbors()
        oxygen_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
        if len(oxygen_neighbors) != 1:
            return False, "Terminal carbon does not have an -OH group"
        # Check that terminal carbon is connected to two hydrogens
        if term_atom.GetTotalNumHs() != 2:
            return False, "Terminal carbon does not have two hydrogens"

    # Check internal carbons
    for atom in chain_atoms[1:-1]:
        # Internal carbons should be connected to one -OH and one hydrogen
        neighbors = atom.GetNeighbors()
        oxygen_neighbors = [nbr for nbr in neighbors if nbr.GetAtomicNum() == 8]
        if len(oxygen_neighbors) != 1:
            return False, "Internal carbon does not have exactly one -OH group"
        # Check that internal carbon is connected to one hydrogen
        if atom.GetTotalNumHs() != 1:
            return False, "Internal carbon does not have one hydrogen"

    # Ensure there are no other functional groups
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in [6, 8, 1]:  # C, O, H
            return False, "Molecule contains atoms other than C, O, and H"
        if atom.GetAtomicNum() == 8:
            # Oxygen atoms should be attached to hydrogen (hydroxyl groups)
            neighbors = atom.GetNeighbors()
            if len(neighbors) != 1 or neighbors[0].GetAtomicNum() != 1:
                return False, "Oxygen atom is not part of a hydroxyl group"

    return True, "Molecule is an alditol with the general formula HOCH2[CH(OH)]nCH2OH"