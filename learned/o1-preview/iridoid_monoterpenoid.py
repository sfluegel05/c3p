"""
Classifies: CHEBI:50563 iridoid monoterpenoid
"""
from rdkit import Chem

def is_iridoid_monoterpenoid(smiles: str):
    """
    Determines if a molecule is an iridoid monoterpenoid based on its SMILES string.

    Iridoids usually consist of a cyclopentane ring fused to a six-membered oxygen heterocycle.
    Secoiridoids are formed by cleavage of a bond in the cyclopentane ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an iridoid monoterpenoid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    if not atom_rings:
        return False, "No rings detected in molecule"

    # Find all five-membered rings
    five_membered_rings = [set(ring) for ring in atom_rings if len(ring) == 5]
    # Find all six-membered rings containing oxygen
    six_membered_oxygen_rings = []
    for ring in atom_rings:
        if len(ring) == 6:
            ring_set = set(ring)
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            atom_symbols = [atom.GetSymbol() for atom in ring_atoms]
            if 'O' in atom_symbols:
                six_membered_oxygen_rings.append(ring_set)

    # Check for fused rings (shared atoms between five- and six-membered rings)
    for five_ring in five_membered_rings:
        for six_ring in six_membered_oxygen_rings:
            shared_atoms = five_ring & six_ring
            if len(shared_atoms) >= 2:
                return True, "Molecule contains cyclopentane ring fused to six-membered oxygen heterocycle"

    # Check for secoiridoid: presence of six-membered oxygen-containing ring
    if six_membered_oxygen_rings:
        return True, "Molecule contains six-membered oxygen heterocycle (possible secoiridoid)"

    return False, "Molecule does not match iridoid monoterpenoid structural features"