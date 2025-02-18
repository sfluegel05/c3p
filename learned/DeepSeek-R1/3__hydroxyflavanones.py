"""
Classifies: CHEBI:48024 3'-hydroxyflavanones
"""
"""
Classifies: 3'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_3__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 3'-hydroxyflavanone based on its SMILES string.
    A 3'-hydroxyflavanone has a flavanone skeleton with a hydroxyl group at the 3' position of the B ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Define flavanone core pattern (dihydrofuran with ketone and phenyl group)
    flavanone_core = Chem.MolFromSmarts("[O]=C1C[C@H](c2ccccc2)CC1")
    if not flavanone_core:
        return False, "Failed to parse flavanone core SMARTS"

    # Check for flavanone structure
    matches = mol.GetSubstructMatches(flavanone_core)
    if not matches:
        return False, "Flavanone core not found"

    # Get the attachment atom (where B ring connects to C ring)
    # In the SMARTS, the phenyl group is atom 3 in the core (0-based index)
    # The core atoms in the match are ordered as per the SMARTS:
    # [O]=C1C[C@H](c2ccccc2)CC1
    # Atoms: 0 (O), 1 (C=O), 2 (C), 3 (C@H), 4 (c2ccccc2), etc.
    # Wait, the SMARTS may have different atom ordering. Let's check:
    # The SMARTS [O]=C1C[C@H](c2ccccc2)CC1 has atoms:
    # 0: O
    # 1: C (double bond to O)
    # 2: C connected to C1
    # 3: C[@H] connected to c2...
    # So the phenyl group is attached to atom 3.

    # For each match, get the attachment atom (atom 3 in the match)
    for match in matches:
        attachment_atom = match[3]
        # Find the B ring (phenyl group connected to attachment_atom)
        # Get the atom's neighbors
        atom = mol.GetAtomWithIdx(attachment_atom)
        neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        # Find the B ring (should be a benzene ring)
        rings = mol.GetRingInfo()
        b_ring = None
        for ring in rings.AtomRings():
            if len(ring) == 6 and attachment_atom in ring:
                # Check if it's a benzene ring (all carbons)
                if all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring):
                    b_ring = ring
                    break
        if not b_ring:
            continue  # No B ring found in this match

        # Check for hydroxyl at 3' position in B ring
        # Find all hydroxyl oxygens in B ring
        hydroxyls = []
        for a in b_ring:
            atom = mol.GetAtomWithIdx(a)
            if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1:
                hydroxyls.append(a)

        # Check if any hydroxyl is two bonds away from attachment_atom
        for oh in hydroxyls:
            path = rdmolops.GetShortestPath(mol, attachment_atom, oh)
            if len(path) == 3:  # Path length 2 bonds (3 atoms)
                return True, "3'-hydroxy group present on B ring"

    return False, "No 3'-hydroxy group found on B ring"