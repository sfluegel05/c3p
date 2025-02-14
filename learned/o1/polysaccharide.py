"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_polysaccharide(smiles: str):
    """
    Determines if a molecule is a polysaccharide based on its SMILES string.
    A polysaccharide is a biomacromolecule consisting of more than ten monosaccharide residues linked glycosidically.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polysaccharide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    # Identify monosaccharide rings (5 or 6 membered rings containing oxygen)
    monosaccharide_rings = []
    for ring in rings:
        # Check ring size (5 or 6 atoms)
        if len(ring) in [5, 6]:
            # Get atomic numbers of atoms in the ring
            atom_nums = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring]
            # Check if ring contains exactly one oxygen atom (common in monosaccharides)
            if atom_nums.count(8) == 1:
                monosaccharide_rings.append(set(ring))

    total_units = len(monosaccharide_rings)

    if total_units == 0:
        return False, "No monosaccharide units found"

    # Build a connectivity graph between rings
    ring_connections = {i: set() for i in range(len(monosaccharide_rings))}
    for bond in mol.GetBonds():
        idx1 = bond.GetBeginAtomIdx()
        idx2 = bond.GetEndAtomIdx()
        # Check if bond connects two different rings
        rings_containing_idx1 = [i for i, ring in enumerate(monosaccharide_rings) if idx1 in ring]
        rings_containing_idx2 = [i for i, ring in enumerate(monosaccharide_rings) if idx2 in ring]
        for ring1 in rings_containing_idx1:
            for ring2 in rings_containing_idx2:
                if ring1 != ring2:
                    # Check if bond is an ether linkage (C-O-C)
                    atom1 = mol.GetAtomWithIdx(idx1)
                    atom2 = mol.GetAtomWithIdx(idx2)
                    if (atom1.GetAtomicNum(), atom2.GetAtomicNum()) in [(6,8),(8,6)]:
                        ring_connections[ring1].add(ring2)
                        ring_connections[ring2].add(ring1)

    # Check if rings are connected via glycosidic linkages to form a polysaccharide chain
    visited = set()
    def dfs(ring_idx):
        visited.add(ring_idx)
        for neighbor in ring_connections[ring_idx]:
            if neighbor not in visited:
                dfs(neighbor)

    # Start DFS from the first ring
    dfs(0)

    # Count connected monosaccharide units
    connected_units = len(visited)

    if connected_units > 10:
        return True, f"Contains {connected_units} monosaccharide units linked glycosidically"
    else:
        return False, f"Contains only {connected_units} monosaccharide units, need more than 10 for a polysaccharide"