"""
Classifies: CHEBI:36141 quinone
"""
from rdkit import Chem

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

    # Define SMARTS pattern for a quinone: a ring with two carbonyl groups attached to sp2 carbons
    quinone_pattern = Chem.MolFromSmarts('[$([cR1]=O)]')

    # Find all carbonyl groups attached to ring carbons
    carbonyl_matches = mol.GetSubstructMatches(quinone_pattern)
    ring_carbonyl_carbons = set()
    for match in carbonyl_matches:
        carbon_idx = match[0]
        oxygen_idx = match[1]
        atom = mol.GetAtomWithIdx(carbon_idx)
        if atom.IsInRing():
            ring_carbonyl_carbons.add(carbon_idx)

    # Check if there are at least two carbonyl groups attached to ring carbons
    if len(ring_carbonyl_carbons) < 2:
        return False, "Less than two carbonyl groups attached to ring carbons"

    # Get all rings in the molecule
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Check for rings that contain at least two carbonyl carbons
    for ring in atom_rings:
        ring_set = set(ring)
        if len(ring_set.intersection(ring_carbonyl_carbons)) >= 2:
            # Check if the ring is fully conjugated (all ring atoms are sp2 hybridized)
            is_conjugated = True
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetHybridization() != Chem.HybridizationType.SP2:
                    is_conjugated = False
                    break
            if is_conjugated:
                return True, "Molecule contains a fully conjugated cyclic diketone characteristic of quinones"

    return False, "Molecule does not contain the characteristic quinone structure"