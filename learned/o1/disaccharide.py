"""
Classifies: CHEBI:36233 disaccharide
"""
from rdkit import Chem

def is_disaccharide(smiles: str):
    """
    Determines if a molecule is a disaccharide based on its SMILES string.
    A disaccharide is a compound in which two monosaccharides are joined by a glycosidic bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a disaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Initialize list to hold sugar rings
    sugar_rings = []

    # Get ring information
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()

    # Map atom index to ring index
    ring_atom_map = {}

    for ring_idx, ring in enumerate(atom_rings):
        ring_size = len(ring)
        if ring_size == 5 or ring_size == 6:
            ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            num_oxygen = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 8)
            num_carbon = sum(1 for atom in ring_atoms if atom.GetAtomicNum() == 6)
            if ring_size == 5 and num_oxygen == 1 and num_carbon == 4:
                sugar_rings.append(set(ring))
                for idx in ring:
                    ring_atom_map[idx] = len(sugar_rings) - 1
            elif ring_size == 6 and num_oxygen == 1 and num_carbon == 5:
                sugar_rings.append(set(ring))
                for idx in ring:
                    ring_atom_map[idx] = len(sugar_rings) - 1

    if len(sugar_rings) != 2:
        return False, f"Found {len(sugar_rings)} sugar rings, need exactly 2"

    # Find glycosidic bonds connecting the two sugar rings
    connected = False

    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()

        if idx1 in ring_atom_map and idx2 in ring_atom_map:
            ring_idx1 = ring_atom_map[idx1]
            ring_idx2 = ring_atom_map[idx2]
            if ring_idx1 != ring_idx2:
                # Atoms in different sugar rings
                atom1 = mol.GetAtomWithIdx(idx1)
                atom2 = mol.GetAtomWithIdx(idx2)
                if atom1.GetAtomicNum() == 8 or atom2.GetAtomicNum() == 8:
                    # One of the atoms is oxygen
                    connected = True
                    break

    if not connected:
        return False, "No glycosidic bond connecting the sugar rings found"

    return True, "Contains two monosaccharide units connected via a glycosidic bond"