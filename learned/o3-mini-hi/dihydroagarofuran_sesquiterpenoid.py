"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: Dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.
The dihydroagarofuran skeleton is expected to be a fully saturated (non‐aromatic) tricyclic core built from a 5–6–5 ring system
that, after removal of decorating substituents (using the Murcko scaffold), contains roughly 15 carbon atoms (allowing 14–16, or one less if one oxygen is present)
and at most one oxygen.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from itertools import combinations

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    Approach:
      1. Parse SMILES.
      2. Compute the Murcko scaffold of the molecule to remove peripheral substituents.
      3. Obtain ring information from the scaffold.
      4. Filter rings to only those that are 5- or 6-membered and non‐aromatic.
      5. Consider all combinations of 3 candidate rings:
           - They must be connected (fused): at least 2 ring pairs share at least 2 atoms.
           - When sorted by ring size they should equal [5, 5, 6].
      6. For a fused candidate triplet, take the union of the ring atom indices (from the scaffold),
         then restrict to atoms that are either carbon or oxygen.
         Finally, count the number of carbons and oxygens.
         The union should have roughly 15 carbons (14–16 if no oxygen is present, or 13–15 if one oxygen is present)
         and at most one oxygen.
    
    Args:
       smiles (str): SMILES string of the molecule.
    
    Returns:
       bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise.
       str: Explanation for the classification decision.
    """
    # Parse the input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Compute the Murcko scaffold to remove most substituents and reveal the core.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not compute Murcko scaffold."
    
    ring_info = scaffold.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in the Murcko scaffold."
    
    # Filter candidate rings: only keep rings of size 5 or 6 that are entirely non‐aromatic
    candidate_rings = []
    for ring in rings:
        if len(ring) not in (5, 6):
            continue
        if any(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        candidate_rings.append(ring)
    
    if len(candidate_rings) < 3:
        return False, (f"Only {len(candidate_rings)} candidate ring(s) (of size 5 or 6 and non‐aromatic) were detected in the scaffold; "
                       "at least 3 are needed for a fused 5–6–5 core.")
    
    # Try each combination of 3 candidate rings.
    for triplet in combinations(candidate_rings, 3):
        # Check connectivity: for fused rings, at least 2 pairs must share at least 2 atoms.
        fused_edges = 0
        for ring_a, ring_b in combinations(triplet, 2):
            if len(set(ring_a).intersection(ring_b)) >= 2:
                fused_edges += 1
        if fused_edges < 2:
            continue

        # Check that the three rings have the expected sizes (order independent): two rings of size 5 and one of size 6.
        ring_sizes = sorted([len(r) for r in triplet])
        if ring_sizes != [5, 5, 6]:
            continue

        # Get the union of atom indices in these rings
        union_atoms = set()
        for ring in triplet:
            union_atoms.update(ring)
        
        # Restrict attention to atoms that are either carbon (6) or oxygen (8)
        core_atoms = []
        for idx in union_atoms:
            atom = scaffold.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() in (6, 8):
                core_atoms.append(atom)
            else:
                core_atoms = None
                break
        if core_atoms is None:
            continue
        
        # Count carbons and oxygens in the candidate core.
        c_count = sum(1 for atom in core_atoms if atom.GetAtomicNum() == 6)
        o_count = sum(1 for atom in core_atoms if atom.GetAtomicNum() == 8)
        
        # Allow at most one oxygen.
        if o_count > 1:
            continue
        # Check that the number of carbon atoms is roughly 15.
        # If one oxygen is present, we allow one less carbon.
        if o_count == 0:
            if not (14 <= c_count <= 16):
                continue
        elif o_count == 1:
            if not (13 <= c_count <= 15):
                continue

        # If we have found a candidate fused 5–6–5 core in the scaffold matching our criteria, return success.
        return True, ("Molecule's scaffold contains a fused 5–6–5 tricyclic non‐aromatic core with "
                      f"{c_count} carbon(s) and {o_count} oxygen(s), consistent with a dihydroagarofuran skeleton.")
    
    # If no candidate core was found, return failure.
    return False, "No fused 5–6–5 tricyclic core (non‐aromatic, with ~15 carbons and at most 1 oxygen) was detected in the molecule's scaffold."