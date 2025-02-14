"""
Classifies: CHEBI:50523 butenolide
"""
"""
Classifies: Butenolide
"""
from rdkit import Chem

def is_butenolide(smiles: str):
    """
    Determines if a molecule is a butenolide based on its SMILES string.
    A butenolide is a gamma-lactone that consists of a 2-furanone skeleton and its substituted derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a butenolide, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    if not ring_info.IsInitialized():
        return False, "Molecule does not contain any rings"

    # Iterate over rings
    for ring_atoms in ring_info.AtomRings():
        # Check if ring is five-membered
        if len(ring_atoms) != 5:
            continue

        ring = [mol.GetAtomWithIdx(idx) for idx in ring_atoms]
        atom_symbols = [atom.GetSymbol() for atom in ring]

        # Check if ring contains exactly one oxygen atom
        if atom_symbols.count('O') != 1:
            continue

        # Check for lactone functionality: oxygen atom in ring is connected to a carbonyl group
        ring_oxygen = None
        for atom in ring:
            if atom.GetSymbol() == 'O':
                ring_oxygen = atom
                break

        if ring_oxygen is None:
            continue  # No oxygen atom found in the ring, should not happen here

        # Check if ring oxygen is connected to a carbon (in ring)
        connected_carbons = [nbr for nbr in ring_oxygen.GetNeighbors() if nbr.GetSymbol() == 'C' and nbr.IsInRing()]
        if not connected_carbons:
            continue  # No in-ring carbon connected to ring oxygen

        ring_carbon = connected_carbons[0]

        # Check if ring carbon is double-bonded to another oxygen (carbonyl group)
        has_carbonyl = False
        for neighbor in ring_carbon.GetNeighbors():
            if neighbor.GetIdx() == ring_oxygen.GetIdx():
                continue  # Skip ring oxygen
            bond = mol.GetBondBetweenAtoms(ring_carbon.GetIdx(), neighbor.GetIdx())
            if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetSymbol() == 'O':
                has_carbonyl = True
                break

        if not has_carbonyl:
            continue  # No carbonyl group found

        # Check for at least one double bond within the ring
        has_double_bond_in_ring = False
        bonds = [mol.GetBondBetweenAtoms(ring_atoms[i], ring_atoms[(i+1)%5]) for i in range(5)]
        for bond in bonds:
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                has_double_bond_in_ring = True
                break

        if not has_double_bond_in_ring:
            continue  # No double bond within the ring

        # All criteria met
        return True, "Contains butenolide ring (gamma-lactone with double bond)"

    return False, "Does not contain butenolide ring"