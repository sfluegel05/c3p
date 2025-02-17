"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: Dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.
The dihydroagarofuran skeleton is expected to be a fully saturated (non‐aromatic) tricyclic core built from a 5–6–5 ring system.
The fused union (after removal of decorating substituents) should contain roughly 15 carbon atoms (allowing 14–16)
and at most 1 oxygen.
"""

from rdkit import Chem
from itertools import combinations

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    
    Process:
      1. Parse the SMILES.
      2. Obtain the ring information.
      3. From all rings, select candidate rings: only those that have 5 or 6 atoms and are non‐aromatic.
         (We no longer restrict by atomic number, to allow the ring to include one oxygen if needed.)
      4. Consider all combinations of 3 candidate rings.
         For a valid fused core, the rings must be connected (at least 2 pairs of rings share at least 2 atoms).
         In addition, the three rings should have sizes that (when sorted) equal [5,5,6].
      5. For each candidate fused triplet, take the union of the ring atom indices. 
         Then, restrict the union to atoms whose atomic numbers are only carbon (6) or oxygen (8).
         Finally, require that the number of carbon atoms is roughly 15 (14–16 allowed) and that there is at most one oxygen.
         
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise.
        str: Explanation for the decision.
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    ring_info = mol.GetRingInfo()
    all_rings = ring_info.AtomRings()
    if not all_rings:
        return False, "No rings detected in the molecule."
    
    # Filter candidate rings:
    # Only rings with 5 or 6 members that are entirely non‐aromatic are kept.
    candidate_rings = []
    for ring in all_rings:
        if len(ring) not in (5, 6):
            continue
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        candidate_rings.append(ring)
        
    if len(candidate_rings) < 3:
        return False, (f"Only {len(candidate_rings)} candidate rings (of size 5 or 6 and non‐aromatic) were detected; "
                       "at least 3 are needed for a fused 5–6–5 system.")
    
    # Try each combination of 3 candidate rings.
    for triplet in combinations(candidate_rings, 3):
        # Check connectivity between rings: for fused rings, at least 2 pairs must share 2 or more atoms.
        fused_edges = 0
        for ring_a, ring_b in combinations(triplet, 2):
            if len(set(ring_a).intersection(ring_b)) >= 2:
                fused_edges += 1
        if fused_edges < 2:
            continue

        # Check that the three rings have the expected sizes (order independent): two rings must be of size 5 and one of size 6.
        ring_sizes = sorted([len(r) for r in triplet])
        if ring_sizes != [5, 5, 6]:
            continue

        # Get the union of atom indices in these rings.
        union_atoms = set()
        for ring in triplet:
            union_atoms.update(ring)
        
        # Now, restrict attention to atoms that are part of the expected core (only C and O allowed)
        core_atoms = []
        for idx in union_atoms:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() in (6, 8):
                core_atoms.append(atom)
            else:
                # If any ring atom is not C or O then skip this candidate.
                core_atoms = None
                break
        if core_atoms is None:
            continue

        # Count carbons and oxygens in the candidate core.
        c_count = sum(1 for atom in core_atoms if atom.GetAtomicNum() == 6)
        o_count = sum(1 for atom in core_atoms if atom.GetAtomicNum() == 8)

        # For a dihydroagarofuran skeleton, the union is expected to have roughly 15 carbons (14-16 allowed).
        # If the union includes one oxygen, then one less carbon is allowed.
        if not (14 <= c_count <= 16):
            continue
        if o_count > 1:
            continue

        # If reached here, we have found a fused 5–6–5 core matching our criteria.
        return True, ("Molecule contains a fused 5–6–5 tricyclic non‐aromatic core with "
                      f"{c_count} carbons and {o_count} oxygen(s), consistent with a dihydroagarofuran skeleton.")
    
    # If no candidate core was found, return failure.
    return False, "No fused 5–6–5 tricyclic core (non‐aromatic, with ~15 carbons and at most 1 oxygen) was detected in the molecule."