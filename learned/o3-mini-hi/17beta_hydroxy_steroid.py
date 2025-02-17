"""
Classifies: CHEBI:35343 17beta-hydroxy steroid
"""
"""
Classifies: 17β-hydroxy steroid
Definition: A 17-hydroxy steroid in which the hydroxy group at position 17 has a beta-configuration.
Improved approach:
  1. Parse the SMILES.
  2. From the ring information, identify “pure‐carbon” rings (i.e. rings in which every atom is a carbon).
  3. Build a graph of connectivity among these rings (two rings are “fused” if they share at least 2 atoms).
  4. Look for a connected fused ring system that contains at least three six–membered rings and one five–membered ring.
  5. Among the –OH groups (defined as an oxygen in –OH attached to an aliphatic carbon that is not carbonyl),
     check if one is attached to a carbon that is part of one of the pure–carbon five–membered rings.
  6. If that carbon has an explicitly defined chirality tag then we take that as a (likely beta) signature.
Note: This approach is heuristic and may not capture all edge cases.
"""

from rdkit import Chem

def is_17beta_hydroxy_steroid(smiles: str):
    """
    Determines whether a molecule is a 17β-hydroxy steroid based on its SMILES string.
    
    The procedure is:
      (a) Parse SMILES; if no rings then return false.
      (b) Identify rings that are entirely carbocyclic (all atoms are carbon).
      (c) Build a fused ring “graph” among these pure-carbon rings (fused if share >= 2 atoms).
      (d) Find any connected (fused) component that has at least three six-membered rings and one five-membered ring.
      (e) In that fused system, look for a –OH group on an aliphatic (sp3) carbon that is a member of one
          of the five-membered rings. Check if that carbon has defined chirality (a proxy for beta stereochemistry).
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): A tuple where the first element is True if the molecule appears to be a 17β–hydroxy steroid,
                   and False otherwise; the second element explains the reasoning.
                   If SMILES parsing fails, returns (False, "Invalid SMILES string").
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No rings found in molecule"
    
    # Identify pure-carbon rings (only carbon atoms; note that in the core steroid nucleus all atoms are carbons)
    pure_carbon_rings = []
    for ring in atom_rings:
        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 for idx in ring):
            pure_carbon_rings.append(set(ring))
    if not pure_carbon_rings:
        return False, "No pure-carbon rings found (steroid nucleus expected to be carbocyclic)"
    
    # Build a connectivity graph among pure-carbon rings.
    n_rings = len(pure_carbon_rings)
    adjacency = {i: set() for i in range(n_rings)}
    for i in range(n_rings):
        for j in range(i+1, n_rings):
            if len(pure_carbon_rings[i].intersection(pure_carbon_rings[j])) >= 2:
                adjacency[i].add(j)
                adjacency[j].add(i)
    
    # Find connected components (fused systems) among the pure-carbon rings.
    visited = set()
    fused_components = []
    for i in range(n_rings):
        if i not in visited:
            stack = [i]
            comp = set()
            while stack:
                cur = stack.pop()
                if cur in comp:
                    continue
                comp.add(cur)
                visited.add(cur)
                for neigh in adjacency[cur]:
                    if neigh not in comp:
                        stack.append(neigh)
            fused_components.append(comp)
    
    # Look for a fused system with at least three six-membered rings and one five-membered ring.
    steroid_component = None
    for comp in fused_components:
        count6 = 0
        count5 = 0
        for idx in comp:
            ring_size = len(pure_carbon_rings[idx])
            if ring_size == 6:
                count6 += 1
            elif ring_size == 5:
                count5 += 1
        if count6 >= 3 and count5 >= 1:
            steroid_component = comp
            break
    if steroid_component is None:
        return False, "Steroid nucleus (fused carbocyclic rings with 3 six-membered and 1 five-membered rings) not found"
    
    # Get all atom indices in the steroid fused system and record the five-membered rings of the fused system.
    steroid_atoms = set()
    five_member_rings = []
    for idx in steroid_component:
        steroid_atoms.update(pure_carbon_rings[idx])
        if len(pure_carbon_rings[idx]) == 5:
            five_member_rings.append(pure_carbon_rings[idx])
    
    # Define a SMARTS for an aliphatic carbon (sp3, not carbonyl) bound to an -OH.
    # We use [CX4;!$([C]=O)][OX2H] to match a carbon (sp3; not immediately double-bonded to O)
    # attached to an -OH.
    oh_pattern = Chem.MolFromSmarts("[CX4;!$([C]=O)][OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    if not oh_matches:
        return False, "No aliphatic -OH group (non-carbonyl) found in the molecule"
    
    # Check if any matched -OH is on an atom that is part of one of the fused pure-carbon five-membered rings.
    for match in oh_matches:
        carbon_idx = match[0]  # carbon attached to the -OH
        if carbon_idx not in steroid_atoms:
            continue  # candidate -OH not on the fused steroid nucleus
        # Is the candidate carbon in one of the pure-carbon five-membered rings?
        in_five = any(carbon_idx in ring for ring in five_member_rings)
        if in_five:
            candidate_atom = mol.GetAtomWithIdx(carbon_idx)
            if candidate_atom.GetChiralTag() != Chem.rdchem.ChiralType.CHI_UNSPECIFIED:
                # Found an -OH on a chiral (and thus stereochemically defined) carbon in a 5-membered ring.
                reason = ("Steroid nucleus (3 six-membered, 1 five-membered, all carbocyclic) found and an -OH "
                          "group is attached on a chiral (likely beta) carbon within one of the five-membered rings "
                          "consistent with a 17β–OH substituent.")
                return True, reason
            else:
                reason = ("Steroid nucleus (3 six-membered, 1 five-membered, all carbocyclic) found and an -OH "
                          "group is attached in a five-membered ring, but the chiral configuration is not defined "
                          "(thus beta configuration is uncertain).")
                return True, reason

    return False, "No -OH group found on a carbon in a pure-carbon five-membered ring of the steroid nucleus (expected for 17β–OH)"

# Example usage:
if __name__ == "__main__":
    # A known positive case:
    test_smiles = "C[C@]12CC[C@H]3[C@@H](CCc4cc(OC(=O)c5ccccc5)ccc34)[C@@H]1CC[C@@H]2O"
    result, reason = is_17beta_hydroxy_steroid(test_smiles)
    print(result, reason)