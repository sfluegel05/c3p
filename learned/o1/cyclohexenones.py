"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: cyclohexenones
"""

from rdkit import Chem

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is any six-membered alicyclic ketone having one double bond in the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cyclohexenone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for cyclohexenone
    # Six-membered non-aromatic ring with one C=O and one C=C in the ring
    cyclohexenone_pattern = Chem.MolFromSmarts(
        "[R6;$([#6]=O)]"
    )
    if cyclohexenone_pattern is None:
        return False, "Invalid SMARTS pattern"

    # Search for the pattern in the molecule
    matches = mol.GetSubstructMatches(cyclohexenone_pattern)

    if not matches:
        return False, "Does not contain a cyclohexenone ring structure"

    # Further filter matches to ensure ring has exactly two double bonds
    ring_info = mol.GetRingInfo()
    for match in matches:
        atom_indices = match

        # Get the ring atoms associated with the matched atoms
        rings = [ring for ring in ring_info.AtomRings() if set(atom_indices).issubset(set(ring))]

        # Iterate over possible rings containing the match
        for ring in rings:
            if len(ring) != 6:
                continue  # Not a six-membered ring

            is_aromatic = any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)
            if is_aromatic:
                continue  # Skip aromatic rings

            # Collect bonds in the ring
            bond_indices = []
            for i in range(len(ring)):
                idx1 = ring[i]
                idx2 = ring[(i + 1) % len(ring)]
                bond = mol.GetBondBetweenAtoms(idx1, idx2)
                if bond is not None:
                    bond_indices.append(bond.GetIdx())

            # Count double bonds in the ring
            double_bond_count = 0
            for bond_idx in bond_indices:
                bond = mol.GetBondWithIdx(bond_idx)
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    double_bond_count += 1

            # Ensure ring has exactly two double bonds (one C=O and one C=C)
            if double_bond_count != 2:
                continue  # Need exactly two double bonds in the ring

            # Check for ketone group (C=O) within the ring
            ketone_in_ring = False
            for bond_idx in bond_indices:
                bond = mol.GetBondWithIdx(bond_idx)
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    begin_atom = bond.GetBeginAtom()
                    end_atom = bond.GetEndAtom()
                    if (begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or \
                       (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6):
                        ketone_in_ring = True
                        break
            if not ketone_in_ring:
                continue  # No ketone group in ring

            # Check for one C=C bond in the ring
            c_c_double_bond = False
            for bond_idx in bond_indices:
                bond = mol.GetBondWithIdx(bond_idx)
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    begin_atom = bond.GetBeginAtom()
                    end_atom = bond.GetEndAtom()
                    if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                        c_c_double_bond = True
                        break
            if not c_c_double_bond:
                continue  # No C=C bond in ring

            # If all conditions are met
            return True, "Contains a six-membered alicyclic ketone ring with one double bond"

    return False, "Does not contain a cyclohexenone ring structure"