"""
Classifies: CHEBI:18154 polysaccharide
"""
"""
Classifies: polysaccharide
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

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

    # Remove hydrogens for simplicity
    mol = Chem.RemoveHs(mol)

    # Define a SMARTS pattern for monosaccharide units (pyranose and furanose forms)
    # This pattern matches a 5 or 6-membered ring with oxygen atoms and hydroxyl groups characteristic of sugars
    monosaccharide_smarts = """
    [
        [C;H1,H2]1
        [C;H1,H2][O][C;H1,H2]
        [C;H1,H2][C;H1,H2]1
    ],
    [
        [C;H1,H2]1
        [C;H1,H2][C;H1,H2][O]
        [C;H1,H2][C;H1,H2]1
    ]
    """

    # Due to the complexity of defining monosaccharides with SMARTS, we'll use a more generic approach
    # Identify all 5 or 6-membered rings with oxygen atoms

    # Get ring information
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()

    monosaccharide_rings = []
    for ring in rings:
        # Check ring size
        if len(ring) in [5, 6]:
            # Get atoms in the ring
            atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
            # Count oxygen atoms in the ring
            oxygens_in_ring = sum(1 for atom in atoms if atom.GetAtomicNum() == 8)
            # Monosaccharide rings usually have one or more oxygen atoms in the ring
            if oxygens_in_ring >= 1:
                monosaccharide_rings.append(set(ring))

    total_units = len(monosaccharide_rings)

    if total_units == 0:
        return False, "No monosaccharide units found"

    # Build a connectivity map between monosaccharide rings via glycosidic linkages
    # Glycosidic linkage pattern: Anomeric carbon (connected to two oxygens) linked to another sugar unit via oxygen

    # Find anomeric carbons: carbons connected to two oxygens
    anomeric_carbons = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 6:
            oxygens = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 8]
            if len(oxygens) == 2:
                anomeric_carbons.append(atom.GetIdx())

    # Build ring connectivity via glycosidic bonds
    ring_connections = {i: set() for i in range(len(monosaccharide_rings))}
    for idx, ring in enumerate(monosaccharide_rings):
        # Find anomeric carbons in the ring
        ring_anomeric_carbons = ring.intersection(anomeric_carbons)
        for ac_idx in ring_anomeric_carbons:
            atom_ac = mol.GetAtomWithIdx(ac_idx)
            for bond in atom_ac.GetBonds():
                nbr = bond.GetOtherAtomIdx(ac_idx)
                if nbr not in ring:
                    # Check if connected to an oxygen that links to another ring
                    atom_nbr = mol.GetAtomWithIdx(nbr)
                    if atom_nbr.GetAtomicNum() == 8:
                        for bond2 in atom_nbr.GetBonds():
                            nbr2 = bond2.GetOtherAtomIdx(nbr)
                            if nbr2 != ac_idx:
                                # Find if nbr2 is in another monosaccharide ring
                                for idx2, ring2 in enumerate(monosaccharide_rings):
                                    if idx2 != idx and nbr2 in ring2:
                                        ring_connections[idx].add(idx2)
                                        break

    # Count connected monosaccharide units
    visited = set()

    def dfs(ring_idx):
        visited.add(ring_idx)
        for neighbor in ring_connections[ring_idx]:
            if neighbor not in visited:
                dfs(neighbor)

    # Perform DFS for all rings to account for branched structures
    for i in range(len(monosaccharide_rings)):
        if i not in visited:
            dfs(i)

    connected_units = len(visited)

    if connected_units > 10:
        return True, f"Contains {connected_units} monosaccharide units linked glycosidically"
    else:
        return False, f"Contains only {connected_units} monosaccharide units, need more than 10 for a polysaccharide"