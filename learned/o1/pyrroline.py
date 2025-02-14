"""
Classifies: CHEBI:23763 pyrroline
"""
from rdkit import Chem

def is_pyrroline(smiles: str):
    """
    Determines if a molecule contains a pyrroline ring based on its SMILES string.
    A pyrroline is a dihydropyrrole, a five-membered heterocyclic ring containing one nitrogen atom and two double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains a pyrroline ring, False otherwise
        str: Reason for classification
    """

    # Parse SMILES string to RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information from the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Flag to indicate if a pyrroline ring is found
    pyrroline_found = False

    # Iterate over each ring in the molecule
    for ring in atom_rings:
        if len(ring) != 5:
            continue  # Skip if ring is not five-membered

        # Get atoms in the ring
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Count nitrogen atoms in the ring
        num_nitrogens = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 7)
        if num_nitrogens != 1:
            continue  # Skip if there is not exactly one nitrogen atom

        # Get bonds in the ring
        ring_bonds = []
        for i in range(len(ring)):
            atom1 = mol.GetAtomWithIdx(ring[i])
            atom2 = mol.GetAtomWithIdx(ring[(i+1)%len(ring)])
            bond = mol.GetBondBetweenAtoms(atom1.GetIdx(), atom2.GetIdx())
            if bond is not None:
                ring_bonds.append(bond)

        # Count double bonds in the ring
        num_double_bonds = sum(1 for bond in ring_bonds if bond.GetBondType() == Chem.BondType.DOUBLE)
        if num_double_bonds != 2:
            continue  # Skip if there are not exactly two double bonds

        # All criteria met; pyrroline ring found
        pyrroline_found = True
        break

    if pyrroline_found:
        return True, "Contains a pyrroline ring (five-membered ring with one nitrogen and two double bonds)"
    else:
        return False, "No pyrroline ring found (requires five-membered ring with one nitrogen and two double bonds)"