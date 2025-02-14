"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_4__hydroxyflavanones(smiles: str):
    """
    Determines if a molecule is a 4'-hydroxyflavanone based on its SMILES string.
    A 4'-hydroxyflavanone is a flavanone with a hydroxy substituent at the 4' position of the B-ring.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 4'-hydroxyflavanone, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a general SMARTS pattern for the flavanone core with atom mapping
    # Atom mappings:
    # [C:1] - carbonyl carbon
    # [O:2] - carbonyl oxygen
    # [C:3]-[C:4]-[O:5] - heterocyclic ring
    # [c:6] - aromatic carbons in A-ring
    # [c:7] - aromatic carbons in B-ring
    flavanone_core_smarts = """
    [$(*(=O)[#6]):1]-[#6:3]-[#6:4]-[#8:5]-1
    -[#6]-[#6]-[#6]-[#6]-[#6]-1
    -[#6]-2:[#6]:[#6]:[#6]:[#6]:[#6]:2
    """

    flavanone_core = Chem.MolFromSmarts(flavanone_core_smarts)
    if flavanone_core is None:
        return None, "Error in flavanone core SMARTS pattern"

    # Find matches for the flavanone core
    flavanone_matches = mol.GetSubstructMatches(flavanone_core, useChirality=False)
    if not flavanone_matches:
        return False, "Flavanone core not found"

    for match in flavanone_matches:
        match_atoms = {}
        for idx, atom_idx in enumerate(match):
            match_atoms[idx+1] = atom_idx  # Atom mappings start from 1

        # Get the atom index for the attachment point of the B-ring
        b_ring_attachment_idx = match_atoms[5]  # Atom [O:5]

        # Find the B-ring (phenyl ring attached to heterocycle)
        aromatic_rings = mol.GetRingInfo().AtomRings()
        b_ring = None
        for ring in aromatic_rings:
            if b_ring_attachment_idx in ring:
                b_ring = ring
                break

        if b_ring is None:
            continue  # B-ring not found in this match

        # Identify the 4' position (para to attachment point)
        # In a phenyl ring, the atom opposite the attachment point is at a path length of 3
        attachment_atom = mol.GetAtomWithIdx(b_ring_attachment_idx)
        for atom_idx in b_ring:
            if atom_idx == b_ring_attachment_idx:
                continue
            path = Chem.rdmolops.GetShortestPath(mol, b_ring_attachment_idx, atom_idx)
            if len(path) == 4:
                # This is the 4' position
                four_prime_atom = mol.GetAtomWithIdx(atom_idx)
                # Check if it has a hydroxy group attached
                has_hydroxy = False
                for neighbor in four_prime_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        has_hydroxy = True
                        break
                if has_hydroxy:
                    return True, "4'-hydroxy group found on B-ring"
        # If no hydroxy group at 4' position in this match, continue to next match
    # No matches with 4'-hydroxy group found
    return False, "4'-hydroxy group on B-ring not found"