"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: Dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.
The dihydroagarofuran skeleton is expected to be a fully saturated (non‐aromatic) tricyclic core 
built from a 5–6–5 ring system whose fused union (after removal of decorating substituents) 
contains roughly 15 carbon atoms (with a margin of ±1) and at most 1 oxygen.
"""

from rdkit import Chem
from itertools import combinations

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    Here we directly use the molecule’s ring information and search for a 5–6–5 fused tricyclic core.
    
    Process:
      1. Parse the SMILES string.
      2. Get the molecule’s ring information.
      3. From all rings, filter for candidate rings that are of size 5 or 6 and whose atoms are non‐aromatic.
         (Also restrict elements to C and O as these are expected in the core.)
      4. From candidate rings, consider all combinations of 3 rings.
         For a valid fused core the three rings must be connected (at least 2 connecting bonds overall)
         and the ring sizes (in some order) should be [5,5,6].
      5. For each candidate fused triplet, take the union of atom indices and count the number
         of carbon and oxygen atoms. The total carbon count should be roughly 15 (14–16 allowed)
         and minimal oxygen (at most 1) should be present.
         
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise.
        str: Explanation for the decision.
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings detected in the molecule."
    
    # Filter candidate rings: only consider 5 or 6 membered rings;
    # also require that none of the atoms is aromatic and that atoms are either C (6) or O (8)
    candidate_rings = []
    for ring in all_rings:
        if len(ring) not in (5,6):
            continue
        valid_ring = True
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            # The dihydroagarofuran core is expected to be fully saturated (non‐aromatic)
            if atom.GetIsAromatic():
                valid_ring = False
                break
            # Allow only Carbon and Oxygen atoms in the core (other atoms likely come from substituents)
            if atom.GetAtomicNum() not in (6,8):
                valid_ring = False
                break
        if valid_ring:
            candidate_rings.append(ring)
    
    if len(candidate_rings) < 3:
        return False, f"Only {len(candidate_rings)} candidate rings (size 5 or 6, non‐aromatic C/O only) were detected; at least 3 are needed in a fused 5–6–5 system."
    
    # Try each combination of 3 candidate rings.
    # For fused rings, an edge exists between two rings if they share at least 2 atoms (common bond).
    valid_core_found = False
    for triplet in combinations(candidate_rings, 3):
        # Check connectivity: in a fused system, at least two pairs of rings must share 2 or more atoms.
        fused_edges = 0
        # Save ring sizes for later check
        ring_sizes = [len(r) for r in triplet]
        
        # For each pair:
        pairs = list(combinations(triplet, 2))
        for ring_a, ring_b in pairs:
            intersection = set(ring_a).intersection(set(ring_b))
            if len(intersection) >= 2:
                fused_edges += 1
        
        # In a 5–6–5 fused system (where the central ring fuses to both outer rings
        # but the outer rings need not be directly fused), we expect at least 2 fused edges.
        if fused_edges < 2:
            continue
        
        # Check that the three rings have the expected sizes:
        if sorted(ring_sizes) != [5,5,6]:
            continue
        
        # Get the union of atom indices for the fused system.
        union_atoms = set()
        for ring in triplet:
            union_atoms.update(ring)
        
        # Count carbons and oxygens in this union.
        c_count = 0
        o_count = 0
        for idx in union_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() == 6:
                c_count += 1
            elif atom.GetAtomicNum() == 8:
                o_count += 1
        
        # Check the core atom counts: allow a small margin (14 to 16 carbons) and at most 1 oxygen.
        if not (14 <= c_count <= 16):
            continue
        if o_count > 1:
            continue
        
        # Lastly, verify that none of these atoms is aromatic.
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in union_atoms):
            continue
        
        # When reaching here the triplet qualifies.
        valid_core_found = True
        return True, ("Molecule contains a fused 5–6–5 tricyclic non‐aromatic core with "
                      f"{c_count} carbons and {o_count} oxygen(s) consistent with a dihydroagarofuran skeleton.")
    
    # If no triplet was accepted, return failure with a cumulative reason.
    return False, "No fused 5–6–5 tricyclic core (non‐aromatic, with ~15 carbons and minimal oxygens) was detected in the molecule."