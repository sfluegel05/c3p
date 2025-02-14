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
            oxygen_atoms = [atom for atom in ring_atoms if atom.GetAtomicNum() == 8]
            # Continue only if there are exactly two oxygen atoms in the ring
            if len(oxygen_atoms) == 2:
                # Find the carbonyl oxygen (double-bonded to carbon)
                carbonyl_oxygen = None
                carbonyl_carbon = None
                ester_oxygen = None
                # Loop over the atoms in the ring to find the carbonyl group
                for atom in ring_atoms:
                    if atom.GetAtomicNum() == 8:
                        for neighbor in atom.GetNeighbors():
                            bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                            if bond.GetBondType() == Chem.BondType.DOUBLE and neighbor.GetAtomicNum() == 6 and neighbor in ring_atoms:
                                carbonyl_oxygen = atom
                                carbonyl_carbon = neighbor
                    # Break if carbonyl oxygen and carbon found
                    if carbonyl_oxygen:
                        break
                # If carbonyl group found, find the ester oxygen
                if carbonyl_oxygen:
                    for atom in oxygen_atoms:
                        if atom != carbonyl_oxygen:
                            ester_oxygen = atom
                            break
                    # Check if ester oxygen is single-bonded to carbons in the ring
                    if ester_oxygen:
                        single_bonds = [bond for bond in ester_oxygen.GetBonds() if bond.GetBondType() == Chem.BondType.SINGLE]
                        ring_bonds = [bond for bond in single_bonds if bond.GetOtherAtom(ester_oxygen) in ring_atoms]
                        # Ester oxygen should have at least two single bonds in the ring
                        if len(ring_bonds) >= 2:
                            gamma_lactone_found = True
                            break

    if gamma_lactone_found:
        return True, "Contains gamma-lactone ring (five-membered lactone ring)"
    else:
        return False, "No gamma-lactone ring found"