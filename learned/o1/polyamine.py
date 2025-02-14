"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: CHEBI:18154 polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic amino compound that contains two or more amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize count of amino groups
    amino_count = 0

    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:  # Only consider nitrogen atoms
            continue

        is_amino = True  # Assume it's an amino group until proven otherwise

        # Exclude nitrogen atoms connected to oxygen or sulfur (e.g., nitro groups)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 8 or nbr.GetAtomicNum() == 16:
                is_amino = False
                break

        # Exclude nitrogen atoms with double or triple bonds
        if any(bond.GetBondType() != Chem.rdchem.BondType.SINGLE for bond in atom.GetBonds()):
            is_amino = False

        # Exclude nitrogen atoms that are part of amides (connected to carbonyl carbon)
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 6:  # Neighbor is carbon
                # Check if carbon is double bonded to oxygen
                for nbr_bond in nbr.GetBonds():
                    if nbr_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        nbr_atom = nbr_bond.GetOtherAtom(nbr)
                        if nbr_atom.GetAtomicNum() == 8:  # Oxygen
                            is_amino = False
                            break
                if not is_amino:
                    break  # Exit early if amide nitrogen

        if is_amino:
            amino_count += 1

    if amino_count >= 2:
        return True, f"Molecule contains {amino_count} amino groups"
    else:
        return False, f"Molecule contains {amino_count} amino group(s), need at least 2"