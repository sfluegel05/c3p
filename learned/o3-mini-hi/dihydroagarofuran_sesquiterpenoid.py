"""
Classifies: CHEBI:71548 dihydroagarofuran sesquiterpenoid
"""
"""
Classifies: Dihydroagarofuran sesquiterpenoid
Definition: Any sesquiterpenoid with a dihydroagarofuran skeleton.
The dihydroagarofuran skeleton is expected to be a fully saturated (non‐aromatic) tricyclic core 
built from a 5–6–5 ring system whose union (after removal of decorating substituents) has roughly 15 carbon atoms 
and very few (at most one) oxygen(s) in the core.
"""

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def is_dihydroagarofuran_sesquiterpenoid(smiles: str):
    """
    Determines if a molecule is a dihydroagarofuran sesquiterpenoid based on its SMILES string.
    The method first extracts the Bemis–Murcko scaffold in order to remove peripheral substituents.
    Then it examines the ring system looking for a fused tricyclic core built (heuristically) from a 5–6–5 system 
    (i.e. with two 5-membered rings and one 6-membered ring fused together). In the fused core, the total number 
    of carbon atoms (ignoring extra oxygens) should be around 15 (we allow a ±1 margin) and there should be at most 
    one oxygen. The backbone should be fully saturated (non‐aromatic).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a dihydroagarofuran sesquiterpenoid, False otherwise
        str: Reason for classification or non-classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
        
    # Extract the Bemis–Murcko scaffold to remove many peripheral substituents.
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None:
        return False, "Could not extract scaffold from molecule."

    # Get ring information from the scaffold
    ring_info = scaffold.GetRingInfo()
    rings = ring_info.AtomRings()
    if not rings:
        return False, "No rings detected in scaffold."
        
    # For a dihydroagarofuran core we expect three fused rings.
    # To narrow our search, first filter rings having size 5 or 6
    candidate_rings = [ring for ring in rings if len(ring) in (5, 6)]
    if len(candidate_rings) < 3:
        return False, f"Only {len(candidate_rings)} rings of size 5 or 6 were found in scaffold; expected 3 for a 5–6–5 tricyclic core."
        
    # Build a connectivity graph for candidate rings.
    # Two rings are considered fused if they share at least 2 atoms (i.e. a common bond).
    ring_count = len(candidate_rings)
    fusion_graph = {i: set() for i in range(ring_count)}
    for i in range(ring_count):
        set_i = set(candidate_rings[i])
        for j in range(i+1, ring_count):
            set_j = set(candidate_rings[j])
            if len(set_i.intersection(set_j)) >= 2:
                fusion_graph[i].add(j)
                fusion_graph[j].add(i)
                
    # Find connected components in the ring fusion graph.
    visited = set()
    components = []
    for i in range(ring_count):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                node = stack.pop()
                if node not in visited:
                    visited.add(node)
                    comp.add(node)
                    stack.extend(fusion_graph[node] - visited)
            components.append(comp)
            
    # Look for a connected component with exactly 3 rings (the candidate dihydroagarofuran core)
    for comp in components:
        if len(comp) == 3:
            # Get the union of atom indices from these three rings
            atom_indices = set()
            ring_sizes = []
            for idx in comp:
                ring = candidate_rings[idx]
                atom_indices.update(ring)
                ring_sizes.append(len(ring))
            # Count carbons and oxygens in the union (from the scaffold)
            c_count = 0
            o_count = 0
            for idx in atom_indices:
                atom = scaffold.GetAtomWithIdx(idx)
                atomic_num = atom.GetAtomicNum()
                if atomic_num == 6:
                    c_count += 1
                elif atomic_num == 8:
                    o_count += 1
            # Heuristic: in an ideal 5–6–5 fused system the union typically contains 15 carbons
            # (we allow a margin of ±1)
            if not (14 <= c_count <= 16):
                reason = f"Tricyclic candidate core has {c_count} carbon atoms; expected ~15."
                continue
            # Ensure that the core is largely hydrocarbonic (allow at most one oxygen)
            if o_count > 1:
                reason = f"Tricyclic candidate core contains {o_count} oxygen atoms; expected at most 1 for the dihydroagarofuran backbone."
                continue
            # Also check that none of the atoms in the union are aromatic.
            if any(scaffold.GetAtomWithIdx(idx).GetIsAromatic() for idx in atom_indices):
                continue

            # Also verify that the three rings have the expected sizes (i.e. two 5-membered and one 6-membered)
            rings_in_comp = [len(candidate_rings[idx]) for idx in comp]
            rings_in_comp.sort()
            if rings_in_comp != [5, 5, 6]:
                reason = f"Ring sizes in candidate core are {rings_in_comp}; expected [5, 5, 6] for a dihydroagarofuran skeleton."
                continue

            return True, "Molecule contains a fused 5–6–5 tricyclic, non‐aromatic core with ~15 carbons consistent with dihydroagarofuran sesquiterpenoid."
            
    return False, "No fused 5–6–5 tricyclic core with ~15 carbons (and minimal oxygens) was detected in the scaffold."

# Example usage (you can uncomment and test with SMILES strings):
# smiles_example = "[C@@]1(C(O[C@@]2([C@](O)([C@@]34...)"  # Replace with a full SMILES
# result, reason = is_dihydroagarofuran_sesquiterpenoid(smiles_example)
# print(result, reason)