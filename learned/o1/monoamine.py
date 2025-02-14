"""
Classifies: CHEBI:63534 monoamine
"""
from rdkit import Chem

def is_monoamine(smiles: str):
    """
    Determines if a molecule is a monoamine based on its SMILES string.
    A monoamine is an aralkylamine compound which contains one amino group connected to an
    aromatic ring by a one- to three-carbon chain. Monoamines are derived from aromatic amino
    acids like phenylalanine, tyrosine, tryptophan, and the thyroid hormones.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a monoamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify potential amine nitrogen atoms
    monoamine_found = False

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:  # Nitrogen atoms only
            continue
        if atom.GetIsAromatic():
            continue  # Exclude aromatic nitrogens
        if any([bond.GetBondType() == Chem.rdchem.BondType.TRIPLE for bond in atom.GetBonds()]):
            continue  # Exclude nitrile nitrogens

        # Exclude amide nitrogens (N connected to C=O)
        is_amide = False
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for bond in neighbor.GetBonds():
                    if bond.GetOtherAtom(neighbor).GetAtomicNum() == 8:  # Oxygen
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetBeginAtomIdx() == neighbor.GetIdx():
                            is_amide = True
                            break
                if is_amide:
                    break
        if is_amide:
            continue

        # Now, find shortest path to any aromatic carbon
        n_idx = atom.GetIdx()
        aromatic_carbons = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 6 and a.GetIsAromatic()]
        if not aromatic_carbons:
            continue  # No aromatic carbons in molecule

        # Compute shortest path lengths to all aromatic carbons
        shortest_paths = []
        for c_idx in aromatic_carbons:
            path = Chem.rdmolops.GetShortestPath(mol, n_idx, c_idx)
            path_length = len(path) - 1  # Number of bonds is one less than number of atoms in path
            shortest_paths.append(path_length)

        min_distance = min(shortest_paths) if shortest_paths else None

        # Check if path length is between 2 and 4 bonds (inclusive)
        if min_distance and 2 <= min_distance <= 4:
            monoamine_found = True
            break  # Stop after finding one match

    if monoamine_found:
        return True, "Contains amino group connected to aromatic ring by a short chain"
    else:
        return False, "Does not contain monoamine functional group"