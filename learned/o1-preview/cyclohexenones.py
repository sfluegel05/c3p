"""
Classifies: CHEBI:48953 cyclohexenones
"""
"""
Classifies: cyclohexenones
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cyclohexenones(smiles: str):
    """
    Determines if a molecule is a cyclohexenone based on its SMILES string.
    A cyclohexenone is defined as any six-membered alicyclic ketone having one double bond in the ring.

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

    # Find all six-membered rings
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    six_membered_rings = [ring for ring in atom_rings if len(ring) == 6]

    if not six_membered_rings:
        return False, "No six-membered rings found"

    # SMARTS pattern for ketone group in ring: [C]=O where C is in ring
    ketone_in_ring_pattern = Chem.MolFromSmarts("[R][C](=O)[R]")

    # Check each six-membered ring
    for ring in six_membered_rings:
        # Check if ring is aliphatic (non-aromatic)
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue  # Skip aromatic rings

        # Get bonds in the ring
        bonds = []
        for i in range(len(ring)):
            bond = mol.GetBondBetweenAtoms(ring[i], ring[(i+1)%len(ring)])
            bonds.append(bond)

        # Count double bonds in the ring
        double_bonds_in_ring = sum(1 for bond in bonds if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE)

        if double_bonds_in_ring != 1:
            continue  # Need exactly one double bond in ring

        # Check for ketone group in the ring
        matches = mol.GetSubstructMatches(ketone_in_ring_pattern)
        for match in matches:
            # Check if both atoms of ketone are in the ring
            if match[0] in ring and match[1] in ring:
                return True, "Contains six-membered alicyclic ring with one double bond and a ketone group in the ring"

    return False, "Does not match criteria for cyclohexenone"

__metadata__ = {
    'chemical_class': {
        'name': 'cyclohexenones',
        'definition': 'Any six-membered alicyclic ketone having one double bond in the ring.'
    }
}