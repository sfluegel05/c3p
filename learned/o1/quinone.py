"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is defined as compounds having a fully conjugated cyclic dione structure,
    derived from aromatic compounds by conversion of an even number of -CH= groups into
    -C(=O)- groups with any necessary rearrangement of double bonds.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Kekulize the molecule to better detect aromaticity
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Chem.KekulizeException:
        # If kekulization fails, proceed without it
        pass

    # Define SMARTS patterns for various quinone types
    quinone_smarts = [
        # Benzoquinone pattern
        'O=C1C=CC=CC1=O',
        # Naphthoquinone pattern
        'O=C1C=CC2=CC=CC=C12',
        # Anthraquinone pattern
        'O=C1C=CC2=CC=CC=C2C1=O',
        # Polycyclic quinone pattern (generalized)
        '[O]=C1[C]=[C][C]=[C][C]=1[O]',
        # Heterocyclic quinone pattern
        '[O]=C1C=CC=NC1=O',
    ]

    # Create pattern molecules
    patterns = [Chem.MolFromSmarts(s) for s in quinone_smarts]

    # Search for quinone patterns
    for pattern in patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule matches quinone SMARTS pattern"

    # Check for fully conjugated cyclic diketone derived from aromatic compounds
    # General approach:
    # - Find rings with two ketone groups
    # - Confirm the ring (or fused ring system) is fully conjugated
    # - Ensure that the ketone positions correspond to the aromatic carbons

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    if not atom_rings:
        return False, "Molecule does not contain any rings"

    # Function to check if a ring system is aromatic before oxidation
    def is_aromatic_before_oxidation(ring_atoms):
        # Create a copy of the molecule
        mol_copy = Chem.RWMol(mol)

        # Reduce the ketones to methylene groups (-C=O to -CH2-)
        for idx in ring_atoms:
            atom = mol_copy.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6 and atom.GetTotalDegree() == 3:
                # Check for C=O group
                for nbr in atom.GetNeighbors():
                    if nbr.GetAtomicNum() == 8 and mol_copy.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx()).GetBondType() == Chem.BondType.DOUBLE:
                        # Replace the double bond O with H (simulating reduction)
                        mol_copy.RemoveBond(atom.GetIdx(), nbr.GetIdx())
                        h_atom = Chem.Atom(1)
                        h_idx = mol_copy.AddAtom(h_atom)
                        mol_copy.AddBond(atom.GetIdx(), h_idx, Chem.BondType.SINGLE)
                        break

        # Check if the ring is aromatic after reduction
        mol_updated = mol_copy.GetMol()
        Chem.SanitizeMol(mol_updated)
        ring_bonds = set()
        for idx in ring_atoms:
            atom = mol_updated.GetAtomWithIdx(idx)
            for bond in atom.GetBonds():
                if bond.GetBeginAtomIdx() in ring_atoms and bond.GetEndAtomIdx() in ring_atoms:
                    ring_bonds.add(bond.GetIdx())
        is_aromatic = True
        for bond_idx in ring_bonds:
            bond = mol_updated.GetBondWithIdx(bond_idx)
            if bond.GetBondType() != Chem.BondType.AROMATIC and bond.GetBondType() != Chem.BondType.DOUBLE:
                is_aromatic = False
                break
        return is_aromatic

    for ring in atom_rings:
        ring_set = set(ring)

        # Find carbonyl carbons in the ring
        carbonyl_carbons = []
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                # Check for C=O group
                is_carbonyl = False
                for nbr in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), nbr.GetIdx())
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_carbonyl = True
                        break
                if is_carbonyl:
                    carbonyl_carbons.append(idx)

        if len(carbonyl_carbons) >= 2 and len(carbonyl_carbons) % 2 == 0:
            # Check if the ring is fully conjugated
            is_conjugated = True
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetHybridization() not in (Chem.HybridizationType.SP2, Chem.HybridizationType.SP):
                    is_conjugated = False
                    break
            if not is_conjugated:
                continue  # Skip this ring

            # Check if the ring was aromatic before oxidation
            if is_aromatic_before_oxidation(ring):
                return True, "Molecule contains a fully conjugated cyclic diketone derived from an aromatic compound"

    return False, "Molecule does not contain the characteristic quinone structure"