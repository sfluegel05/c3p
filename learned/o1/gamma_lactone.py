"""
Classifies: CHEBI:37581 gamma-lactone
"""
"""
Classifies: gamma-lactone
"""

from rdkit import Chem

def is_gamma_lactone(smiles: str):
    """
    Determines if a molecule is a gamma-lactone based on its SMILES string.
    A gamma-lactone is defined as a lactone (cyclic ester) having a five-membered ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a gamma-lactone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    # Get list of rings (list of atom index tuples)
    rings = ri.AtomRings()

    gamma_lactone_found = False

    # Loop over the rings
    for ring in rings:
        if len(ring) == 5:
            # Collect atoms in the ring
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Count number of oxygen atoms in the ring
            n_oxygen_atoms = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            # Continue only if there is exactly one oxygen atom in the ring
            if n_oxygen_atoms == 1:
                # Check for carbonyl group in the ring (C=O)
                for atom in ring_atoms:
                    if atom.GetAtomicNum() == 6:  # Carbon atom
                        # Check if atom has a double bond to oxygen
                        for bond in atom.GetBonds():
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                neighbor = bond.GetOtherAtom(atom)
                                if neighbor.GetAtomicNum() == 8 and neighbor.IsInRing():
                                    # Found carbonyl oxygen in ring
                                    carbonyl_carbon = atom
                                    carbonyl_oxygen = neighbor
                                    # Find the ring oxygen (single bonded oxygen)
                                    ring_oxygen = next((a for a in ring_atoms if a.GetAtomicNum() == 8 and a != carbonyl_oxygen), None)
                                    if ring_oxygen is not None:
                                        # Check if ring oxygen is single bonded to carbonyl carbon (forming ester linkage)
                                        bond_to_oxygen = mol.GetBondBetweenAtoms(carbonyl_carbon.GetIdx(), ring_oxygen.GetIdx())
                                        if bond_to_oxygen and bond_to_oxygen.GetBondType() == Chem.BondType.SINGLE:
                                            # Confirm both atoms are in the same five-membered ring
                                            if ri.IsAtomInRingOfSize(carbonyl_carbon.GetIdx(), 5) and ri.IsAtomInRingOfSize(ring_oxygen.GetIdx(), 5):
                                                gamma_lactone_found = True
                                                break
                            if gamma_lactone_found:
                                break
                    if gamma_lactone_found:
                        break
                if gamma_lactone_found:
                    break

    if gamma_lactone_found:
        return True, "Contains gamma-lactone ring (five-membered lactone ring)"
    else:
        return False, "No gamma-lactone ring found"