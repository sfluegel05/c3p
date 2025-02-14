"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:17164 primary amine
"""

from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is a compound where a nitrogen atom is bonded to one hydrocarbyl group and two hydrogens,
    and the nitrogen is not connected to a carbonyl group (excluding amides).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Nitrogen
            if atom.GetFormalCharge() != 0:
                continue  # Skip charged nitrogen atoms

            if atom.GetDegree() != 1:
                continue  # Nitrogen must be connected to exactly one heavy atom

            if atom.GetTotalNumHs() != 2:
                continue  # Nitrogen must have two hydrogens

            # Get neighbor atom (only one heavy atom neighbor)
            neighbor_atom = atom.GetNeighbors()[0]
            if neighbor_atom.GetAtomicNum() != 6:
                continue  # Neighbor is not carbon

            # Check that the neighbor carbon is not a carbonyl carbon (exclude amides)
            is_amide = False
            for bond in neighbor_atom.GetBonds():
                bonded_atom = bond.GetOtherAtom(neighbor_atom)
                if bonded_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                    is_amide = True
                    break
            if is_amide:
                continue  # Exclude amide nitrogen

            # Check that nitrogen is not part of a nitrile group (C#N)
            if neighbor_atom.GetBondBetweenAtoms(neighbor_atom.GetIdx(), atom.GetIdx()).GetBondType() == Chem.BondType.TRIPLE:
                continue  # Exclude nitriles

            # If all conditions are met, it's a primary amine
            return True, "Primary amine group found"

    return False, "No primary amine group found"