"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is any oxolane (tetrahydrofuran ring with 4 carbons and 1 oxygen)
    having an oxo substituent (C=O) at any position on the tetrahydrofuran ring. This includes lactones
    where the carbonyl is part of the ring, and compounds where a C=O is attached to 
    any atom of the ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES into RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Iterate over all rings in the molecule
    for idx_ring, ring in enumerate(atom_rings):
        # Check if the ring is 5-membered
        if len(ring) != 5:
            continue

        # Initialize counters
        oxygen_count = 0
        carbon_count = 0
        other_atom = False

        # Collect ring atoms
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]

        # Count atoms in the ring
        for atom in ring_atoms:
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                carbon_count += 1
            else:
                # Found an atom other than carbon or oxygen in ring
                other_atom = True
                break
            # Check if atom is aromatic
            if atom.GetIsAromatic():
                other_atom = True
                break
        if other_atom:
            continue

        # Continue only if ring has 1 oxygen, 4 carbons
        if oxygen_count != 1 or carbon_count != 4:
            continue

        # Check for unsaturation within the ring (exclude bonds that are not single)
        bonds_in_ring = [bond for bond in mol.GetBonds() if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring]
        ring_unsaturation = False
        for bond in bonds_in_ring:
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                ring_unsaturation = True
                break
        if ring_unsaturation:
            continue

        # Check that all ring atoms are sp3 hybridized (saturated)
        ring_atoms_sp3 = True
        for atom in ring_atoms:
            if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
                ring_atoms_sp3 = False
                break
        if not ring_atoms_sp3:
            continue

        # Check for presence of C=O bond either in the ring (lactone) or attached to ring atoms
        found_C_O_double = False

        # First, check for C=O bonds in the ring (lactone)
        for bond in bonds_in_ring:
            if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                begin_atom = bond.GetBeginAtom()
                end_atom = bond.GetEndAtom()
                if ((begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 8) or
                    (begin_atom.GetAtomicNum() == 8 and end_atom.GetAtomicNum() == 6)):
                    found_C_O_double = True
                    break

        # If not found in ring, check for C=O attached to ring atoms
        if not found_C_O_double:
            for atom in ring_atoms:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in ring:
                        continue  # Skip ring atoms
                    bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                    if bond is None:
                        continue
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        if neighbor.GetAtomicNum() == 8:
                            # Found exocyclic C=O attached to ring atom
                            found_C_O_double = True
                            break
                if found_C_O_double:
                    break

        if found_C_O_double:
            # Match found
            return True, "Contains tetrahydrofuran ring with oxo substituent at ring position"
        else:
            continue  # No C=O found, check next ring

    # No matching ring found
    return False, "No tetrahydrofuranone ring found"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:[ID]',
        'name': 'tetrahydrofuranone',
        'definition': 'Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.'
    }
}