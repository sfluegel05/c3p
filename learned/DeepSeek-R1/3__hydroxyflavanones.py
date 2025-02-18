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

    # Generalized flavanone core SMARTS without stereochemistry requirements
    # Matches O=C1CC(Ph)CC1 where Ph is any aromatic ring (B ring)
    flavanone_core = Chem.MolFromSmarts("[O]=C1CC(c2ccccc2)CC1")
    if not flavanone_core:
        return False, "Failed to parse flavanone core SMARTS"

    # Check for flavanone structure
    matches = mol.GetSubstructMatches(flavanone_core)
    if not matches:
        return False, "Flavanone core not found"

    # Iterate through all core matches (handle multiple possible conformations)
    for match in matches:
        try:
            # The phenyl group (B ring) is attached to atom index 3 in the SMARTS pattern
            # SMARTS breakdown: [O]=C1(atom0)-C(atom1)-C(atom2)-(c2ccccc2)(atom3)-C(atom4)-C(atom5)-1
            attachment_atom_idx = match[3]  # Atom connecting to B ring
            
            # Get the actual B ring atom connected to the flavanone core
            core_atom = mol.GetAtomWithIdx(attachment_atom_idx)
            b_ring_atoms = [n for n in core_atom.GetNeighbors() if n.GetAtomicNum() == 6]
            
            for b_ring_atom in b_ring_atoms:
                # Verify B ring is a benzene ring
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if b_ring_atom.GetIdx() in ring and len(ring) == 6:
                        # Check if all atoms in ring are carbons (allowing substituents)
                        if all(mol.GetAtomWithIdx(a).GetAtomicNum() == 6 for a in ring):
                            b_ring = ring
                            # Check hydroxyl positions in B ring
                            for atom_idx in b_ring:
                                atom = mol.GetAtomWithIdx(atom_idx)
                                if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() >= 1:
                                    # Calculate path from core attachment to hydroxyl oxygen
                                    path = rdmolops.GetShortestPath(mol, attachment_atom_idx, atom_idx)
                                    # Path length 4 means 3 bonds between core attachment and hydroxyl (positions 1->2->3)
                                    if len(path) == 4:
                                        return True, "3'-hydroxy group present on B ring"
        except IndexError:
            continue  # Handle potential index errors in match processing

    return False, "No 3'-hydroxy group found on B ring"