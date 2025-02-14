"""
Classifies: CHEBI:140331 4'-hydroxyflavanones
"""
from rdkit import Chem

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

    # Define the flavanone core SMARTS pattern with atom mapping
    flavanone_core_smarts = """
    [#6]-1(=[#8])                      # carbonyl group
    -[#6]-[#6]-[#8]-[#6]-1             # heterocyclic ring with oxygen
    -[#6]-2:[#6]:[#6]:[#6]:[#6]:[#6]:2 # A-ring (aromatic ring fused to heterocycle)
    """

    flavanone_core = Chem.MolFromSmarts(flavanone_core_smarts)
    if flavanone_core is None:
        return None, "Error in flavanone core SMARTS pattern"

    # Find matches for the flavanone core
    flavanone_matches = mol.GetSubstructMatches(flavanone_core, useChirality=False)
    if not flavanone_matches:
        return False, "Flavanone core not found"

    for match in flavanone_matches:
        # Map the indices from the SMARTS pattern to the molecule
        # Assuming the first atom in SMARTS is index 0, adjust if needed
        carbonyl_c_idx = match[0]
        oxygen_idx = match[1]
        hetero_ring_c_idx = match[2]
        o_in_ring_idx = match[3]
        fused_ring_c_idx = match[4]

        # Identify the B-ring attached to the heterocycle
        attachment_atom = mol.GetAtomWithIdx(o_in_ring_idx)
        neighbors = [nbr.GetIdx() for nbr in attachment_atom.GetNeighbors() if nbr.GetIdx() != hetero_ring_c_idx]
        if not neighbors:
            continue  # No B-ring attached

        b_ring_start_idx = neighbors[0]
        # Get the B-ring (phenyl ring attached via oxygen)
        rings = mol.GetRingInfo().AtomRings()
        b_ring = None
        for ring in rings:
            if b_ring_start_idx in ring and len(ring) == 6:
                b_ring = ring
                break
        if b_ring is None:
            continue  # B-ring not found

        # Identify the 4' position (para to attachment point)
        attachment_idx = o_in_ring_idx
        b_ring_atoms = set(b_ring)
        for atom_idx in b_ring:
            if atom_idx == b_ring_start_idx:
                continue
            path = Chem.rdmolops.GetShortestPath(mol, b_ring_start_idx, atom_idx)
            if len(path) == 4:
                # This is the 4' position
                four_prime_atom = mol.GetAtomWithIdx(atom_idx)
                # Check if a hydroxy group is attached
                has_hydroxy = False
                for neighbor in four_prime_atom.GetNeighbors():
                    if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() == 1:
                        has_hydroxy = True
                        break
                if has_hydroxy:
                    return True, "4'-hydroxy group found on B-ring"
        # If no hydroxy at 4' position in this match, continue to next match

    # No matches with 4'-hydroxy group found
    return False, "4'-hydroxy group on B-ring not found"