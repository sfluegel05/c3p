"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is any oxolane (tetrahydrofuran ring) having an oxo substituent (C=O) at any position on the ring.

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
    for ring in atom_rings:
        # Check if the ring is 5-membered
        if len(ring) != 5:
            continue

        # Initialize counters
        oxygen_count = 0
        carbon_count = 0
        other_atom = False

        # Count atoms in the ring
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            atomic_num = atom.GetAtomicNum()
            if atomic_num == 8:
                oxygen_count += 1
            elif atomic_num == 6:
                carbon_count += 1
            else:
                # Found an atom other than carbon or oxygen in ring
                other_atom = True
                break

        # Continue only if ring has 1 oxygen, 4 carbons, and no other atoms
        if other_atom or oxygen_count !=1 or carbon_count !=4:
            continue

        # Check for oxo substituent (C=O) on ring carbons
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:  # Carbon atom
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() in ring:
                        continue  # Skip atoms in the ring
                    if neighbor.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx())
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            # Found ring carbon double-bonded to oxygen
                            return True, "Contains tetrahydrofuran ring with oxo substituent"

    # No matching ring found
    return False, "No tetrahydrofuranone ring found"

__metadata__ = {   
    'chemical_class': {   
        'id': 'CHEBI:[ID]',
        'name': 'tetrahydrofuranone',
        'definition': 'Any oxolane having an oxo- substituent at any position on the tetrahydrofuran ring.'
    }
}