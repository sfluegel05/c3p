"""
Classifies: CHEBI:15889 sterol
"""
"""
Classifies a sterol as defined by:
  "Any 3-hydroxy steroid whose skeleton is closely related to cholestan-3-ol 
   (additional carbon atoms may be present in the side chain)."

Our improved approach:
  1. Parse the SMILES string.
  2. Get all rings from the molecule and build a connectivity graph between rings 
     (two rings are connected if they share at least one atom).
  3. Identify the largest fused (connected) ring system. In a sterol nucleus, this should
     be the tetracyclic system (at least 4 fused rings).
  4. Check that the fused ring system is carbon-rich (using a threshold on the number 
     of carbon atoms; we now require at least 16 carbons) and that most of its atoms are carbon.
  5. Look for at least one hydroxyl group (–OH) attached directly to one of the ring carbons.
"""

from rdkit import Chem

def is_sterol(smiles: str):
    """
    Determines if a molecule is a sterol.
    A sterol is defined as any 3‑hydroxy steroid which has a fused (tetracyclic) ring system
    closely related to cholestan-3-ol. Typically, the steroid nucleus is composed of three
    six‑membered rings and one five‑membered ring that are fused, and a hydroxyl group
    must be attached directly to one of the nucleus carbons.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is classified as a sterol, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Get all rings from the molecule
    ring_info = mol.GetRingInfo()
    rings = ring_info.AtomRings()  # list of tuples (each tuple contains atom indices)
    if not rings:
        return False, "No rings found in molecule"
    
    # Build a connectivity graph between rings: two rings are connected if they share at least one atom.
    num_rings = len(rings)
    ring_adj = {i: set() for i in range(num_rings)}
    for i in range(num_rings):
        for j in range(i+1, num_rings):
            # if rings i and j share at least one atom, add an edge
            if set(rings[i]).intersection(rings[j]):
                ring_adj[i].add(j)
                ring_adj[j].add(i)
    
    # Identify connected components among the rings (each component is a fused ring system).
    visited = set()
    components = []
    for i in range(num_rings):
        if i in visited:
            continue
        comp = set()
        stack = [i]
        while stack:
            node = stack.pop()
            if node in visited:
                continue
            visited.add(node)
            comp.add(node)
            stack.extend(ring_adj[node] - visited)
        components.append(comp)
    
    # Choose the connected component that contains the largest number of unique atoms.
    best_atoms = set()
    best_component = None
    for comp in components:
        comp_atoms = set()
        for idx in comp:
            comp_atoms.update(rings[idx])
        if len(comp_atoms) > len(best_atoms):
            best_atoms = comp_atoms
            best_component = comp
    
    # In a typical steroid nucleus, there are 4 fused rings.
    if best_component is None or len(best_component) < 4:
        return False, "No fused tetracyclic (4-ring) system found in molecule"
    
    # Count the number of carbon atoms and total atoms in the largest fused ring system.
    carbon_count = 0
    for atom_idx in best_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        if atom.GetAtomicNum() == 6:
            carbon_count += 1
    total_atoms_in_nucleus = len(best_atoms)
    
    # We relax the threshold slightly (e.g., some steroids like B-norcholesterol may have 16 carbons),
    # and also require that the majority of the nucleus are carbons.
    if carbon_count < 16:
        return False, (f"Fused ring system contains too few carbon atoms (found {carbon_count}, need >= 16)")
    if (carbon_count / total_atoms_in_nucleus) < 0.75:
        return False, "Fused ring system contains too many heteroatoms to be considered a steroid nucleus"
    
    # Look for a hydroxyl group. We use a SMARTS for hydroxyl groups: oxygen with one hydrogen.
    hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)
    
    # Check that at least one hydroxyl group is directly attached to a carbon atom in the fused nucleus.
    found_ring_hydroxyl = False
    for match in hydroxyl_matches:
        o_idx = match[0]  # the oxygen atom of the -OH group
        o_atom = mol.GetAtomWithIdx(o_idx)
        for neighbor in o_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in best_atoms:
                found_ring_hydroxyl = True
                break
        if found_ring_hydroxyl:
            break
    if not found_ring_hydroxyl:
        return False, ("No hydroxyl group found directly attached to the fused tetracyclic ring system")
    
    return True, ("Molecule has a fused tetracyclic steroid nucleus (with most atoms being carbon and at least 16 ring carbons) "
                  "and a hydroxyl group attached to one of its ring carbons; "
                  "thus, it is classified as a sterol.")

# Example usage:
if __name__ == "__main__":
    # Test using one of the example SMILES for a sterol
    test_smiles = "C[C@@H](CO)[C@@H](O)CC[C@@H](C)[C@H]1CC[C@H]2[C@@H]3CC=C4C[C@@H](O)CC[C@]4(C)[C@H]3CC[C@]12C"
    result, reason = is_sterol(test_smiles)
    print("Is sterol:", result)
    print("Reason:", reason)