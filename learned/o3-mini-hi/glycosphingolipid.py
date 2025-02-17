"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid

A glycosphingolipid is defined as a glycolipid that is a carbohydrate‐containing derivative
of a sphingoid or ceramide in which the carbohydrate residue is attached (typically via an O‐glycosidic linkage)
to the O-1 position of the sphingoid.

Improved heuristic strategy:
  1. Identify candidate sugar rings: a 5- or 6-membered ring containing exactly one oxygen.
  2. Identify a sphingoid/ceramide backbone via a relaxed SMARTS pattern.
     Here we use "[#6]C(=O)N[C@@H](CO)" so that the fatty-acid carbonyl and the subsequent amide bearing sphingoid
     alcohol are captured.
  3. For each candidate sugar ring, examine each carbon in the ring (a candidate anomeric carbon) that bears an
     exocyclic oxygen. From that oxygen, perform a breadth-first search (with a maximum depth of 3 bonds) while
     excluding atoms of the sugar ring. If any search touches a sphingoid atom we call it a glycosidic linkage.
  4. Require that the molecule contains at least 20 carbon atoms.
  
If any check fails an explanatory message is returned.
"""

from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       bool: True if the molecule is classified as a glycosphingolipid, False otherwise.
       str: Reason for the classification.
    """
    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check molecule has at least 20 carbon atoms
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons to be a glycosphingolipid"

    # Step 1: Identify candidate sugar rings.
    # Heuristically, a sugar ring is defined as a 5- or 6-membered ring containing exactly one oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # Will store sets of atom indices representing candidate sugar rings.
    for ring in ring_info:
        if len(ring) in (5,6):
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count == 1:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety detected"

    # Step 2: Identify sphingoid/ceramide backbone.
    # Use a refined SMARTS pattern which captures an amide-bearing carbonyl connected to a sphingoid alcohol.
    sphingoid_smarts = "[#6]C(=O)N[C@@H](CO)"
    sphingo_pattern = Chem.MolFromSmarts(sphingoid_smarts)
    sphingo_matches = mol.GetSubstructMatches(sphingo_pattern)
    if not sphingo_matches:
        return False, "No sphingoid/ceramide backbone detected"

    # Collect all atom indices from any sphingoid match.
    sphingo_atoms = set()
    for match in sphingo_matches:
        sphingo_atoms.update(match)
    
    # Helper function: perform breadth-first search (BFS) starting from a given atom (start_atom)
    # up to a given maximum depth (here 3), while excluding atoms in the 'exclude' set.
    # If any atom in the BFS is in sphingo_atoms, return True.
    def search_link(start_atom, exclude, max_depth=3):
        visited = {start_atom.GetIdx()}
        frontier = [(start_atom, 0)]
        while frontier:
            current_atom, depth = frontier.pop(0)
            if depth >= max_depth:
                continue
            for nbr in current_atom.GetNeighbors():
                nid = nbr.GetIdx()
                if nid in exclude or nid in visited:
                    continue
                if nid in sphingo_atoms:
                    return True
                visited.add(nid)
                frontier.append((nbr, depth + 1))
        return False

    # Step 3: Try to find a glycosidic linkage:
    # For each candidate sugar ring, examine its carbon atoms (possible anomeric centers)
    # that have exocyclic oxygen substituents. From each such oxygen, perform BFS (depth up to 3)
    # excluding the sugar ring atoms.
    glyco_link_found = False
    for ring in sugar_rings:
        for idx in ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6:  # Only consider carbons in the ring.
                continue
            # Look at each neighbor of the candidate carbon that is not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Look for exocyclic oxygen substituents (only single bonds are expected in glycosidic bonds)
                if nbr.GetAtomicNum() == 8:
                    # Occasionally the oxygen may directly be bonded to a sphingoid atom.
                    direct_found = any(nbr2.GetIdx() in sphingo_atoms and nbr2.GetIdx() != atom.GetIdx()
                                       for nbr2 in nbr.GetNeighbors())
                    if direct_found:
                        glyco_link_found = True
                        break
                    # Otherwise, perform BFS (up to depth 3) starting from this oxygen.
                    exclude_set = set(ring)
                    exclude_set.add(atom.GetIdx())
                    if search_link(nbr, exclude_set, max_depth=3):
                        glyco_link_found = True
                        break
            if glyco_link_found:
                break
        if glyco_link_found:
            break

    if not glyco_link_found:
        return False, "Sugar and sphingoid moieties are not glycosidically linked (even allowing a bridging group)"
    
    return True, ("Molecule contains a glycosphingolipid structure: a sugar ring is glycosidically linked (directly or via "
                  "a short bridging group) to a sphingoid/ceramide backbone")

# Example usage (for testing purposes):
# smiles_example = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_glycosphingolipid(smiles_example)
# print(result, reason)