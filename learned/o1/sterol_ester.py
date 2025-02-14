"""
Classifies: CHEBI:35915 sterol ester
"""
"""
Classifies: CHEBI:35366 sterol ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sterol_ester(smiles: str):
    """
    Determines if a molecule is a sterol ester based on its SMILES string.
    A sterol ester is a steroid ester obtained by formal condensation of the carboxy group of any carboxylic acid with the 3-hydroxy group of a sterol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sterol ester, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts('[CX3](=O)[OX2H0][#6]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No ester group found"

    # Get ring information
    ri = mol.GetRingInfo()
    atom_rings = ri.AtomRings()

    # Check for rings of size 5 and 6
    ring_sizes = [len(ring) for ring in atom_rings]
    num_6_membered_rings = ring_sizes.count(6)
    num_5_membered_rings = ring_sizes.count(5)

    if num_6_membered_rings < 3 or num_5_membered_rings < 1:
        return False, "Not enough rings of sizes 6 and 5 for steroid nucleus"

    # Build ring adjacency list
    ring_adj_list = {}
    for i in range(len(atom_rings)):
        ring_adj_list[i] = set()
    for i in range(len(atom_rings)):
        for j in range(i+1, len(atom_rings)):
            if set(atom_rings[i]) & set(atom_rings[j]):
                ring_adj_list[i].add(j)
                ring_adj_list[j].add(i)

    # Find fused ring systems
    def find_fused_ring_system(ring_idx, visited=None):
        if visited is None:
            visited = set()
        visited.add(ring_idx)
        for neighbor in ring_adj_list[ring_idx]:
            if neighbor not in visited:
                find_fused_ring_system(neighbor, visited)
        return visited

    for i in range(len(atom_rings)):
        fused_system = find_fused_ring_system(i)
        if len(fused_system) >= 4:
            sizes = [ring_sizes[r] for r in fused_system]
            if sizes.count(6) >= 3 and sizes.count(5) >= 1:
                # Potential steroid nucleus found
                # Now check if ester oxygen is connected to any atom in the fused ring system
                fused_atoms = set()
                for r in fused_system:
                    fused_atoms.update(atom_rings[r])
                for ester_match in ester_matches:
                    ester_oxygen = ester_match[1]
                    # Check if ester oxygen is connected to fused ring system
                    oxygen_atom = mol.GetAtomWithIdx(ester_oxygen)
                    neighbors = oxygen_atom.GetNeighbors()
                    for neighbor in neighbors:
                        if neighbor.GetIdx() in fused_atoms:
                            return True, "Ester group connected to steroid nucleus"
                return False, "Ester group not connected to steroid nucleus"

    return False, "No steroid-like fused ring system found"