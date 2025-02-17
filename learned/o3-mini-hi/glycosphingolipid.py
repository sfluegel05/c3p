"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid

A glycosphingolipid is defined as a glycolipid that is a carbohydrate‐containing derivative
of a sphingoid or ceramide, in which the carbohydrate residue is attached (typically via an O‐glycosidic linkage)
to the O-1 position of the sphingoid.

Our improved heuristic strategy is:
  1. Identify candidate sugar rings using RDKit ring info: a sugar ring is approximated as a 5- or 6-membered ring 
     that contains exactly one oxygen.
  2. Identify a sphingoid/ceramide backbone using a relaxed SMARTS (e.g. "C(=O)NC(CO)") and collect the matching atoms.
  3. Look for a glycosidic linkage between some atom in one of the sugar rings and the sphingoid fragment. 
     For each carbon in a sugar ring (candidate anomeric carbon) with at least one exocyclic oxygen (neighbor not in the ring),
     perform a short breadth-first search (up to depth 2) on that exocyclic oxygen – excluding atoms of the sugar ring – 
     to see if a connection to a sphingoid atom can be found. This allows for one or two “bridging” atoms (such as a phosphorus or additional oxygen).
  4. To reduce false positives, require that the molecule contains at least 20 carbon atoms.
  
If any of the above checks fail, the function returns False with a reason.
"""

from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    
    Args:
       smiles (str): SMILES string of the molecule
       
    Returns:
       bool: True if molecule is classified as a glycosphingolipid, False otherwise.
       str: Reason for the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check that the molecule contains a sufficient number of carbons
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons to be a glycosphingolipid"
    
    # Step 1: Identify candidate sugar rings.
    # Heuristically, a sugar ring is a 5- or 6-membered ring containing exactly one oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # list of sets of atom indices for candidate sugar rings
    for ring in ring_info:
        if len(ring) in (5, 6):
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count == 1:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety detected"
    
    # Step 2: Identify the sphingoid/ceramide backbone.
    # We use a relaxed SMARTS pattern for a ceramide-like fragment.
    sphingoid_smarts = "C(=O)NC(CO)"
    sphingo_pattern = Chem.MolFromSmarts(sphingoid_smarts)
    sphingo_matches = mol.GetSubstructMatches(sphingo_pattern)
    if not sphingo_matches:
        return False, "No sphingoid/ceramide backbone detected"
    # Collect all atom indices from any sphingoid match.
    sphingo_atoms = set()
    for match in sphingo_matches:
        sphingo_atoms.update(match)
    
    # Helper function: perform BFS (up to max_depth) from candidate oxygen (excluding sugar ring atoms given by 'exclude')
    # to determine if any reached atom is in sphingo_atoms.
    def search_link(start_atom, exclude, max_depth=2):
        visited = set([start_atom.GetIdx()])
        frontier = [(start_atom, 0)]
        while frontier:
            curr_atom, depth = frontier.pop(0)
            if depth > max_depth:
                continue
            # Check all neighbors (except those in the excluded set)
            for nbr in curr_atom.GetNeighbors():
                nid = nbr.GetIdx()
                if nid in exclude or nid in visited:
                    continue
                # If we find a sphingoid atom in the search, return True.
                if nid in sphingo_atoms:
                    return True
                visited.add(nid)
                frontier.append((nbr, depth + 1))
        return False
    
    # Step 3: Look for a glycosidic linkage.
    # For each candidate sugar ring, look through its carbon atoms (potential anomeric carbons)
    # that have "exocyclic" oxygen substituents. For each such oxygen, allow a short connection search
    # (up to 2 bonds) toward the sphingoid fragment.
    glyco_link_found = False
    for ring in sugar_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Consider only carbon atoms for anomeric candidates.
            if atom.GetAtomicNum() != 6:
                continue
            # For each neighbor of the candidate carbon that is NOT in the sugar ring (exocyclic),
            # and is an oxygen, search for a connection to the sphingoid part.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue  # skip atoms in the sugar ring
                if nbr.GetAtomicNum() == 8:
                    # Option A: Check direct connection to a sphingoid atom.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() != atom_idx and nbr2.GetIdx() in sphingo_atoms:
                            glyco_link_found = True
                            break
                    if glyco_link_found:
                        break
                    # Option B: Look one or two bonds away via a BFS search.
                    # Exclude the sugar ring (the union of all atoms in this ring plus the current candidate sugar carbon)
                    exclude_set = set(ring)
                    exclude_set.add(atom_idx)
                    if search_link(nbr, exclude_set, max_depth=2):
                        glyco_link_found = True
                        break
            if glyco_link_found:
                break
        if glyco_link_found:
            break

    if not glyco_link_found:
        return False, "Sugar and sphingoid moieties are not glycosidically linked (even allowing a bridging group)"
    
    return True, "Molecule contains a glycosphingolipid structure: a sugar ring is glycosidically linked (directly or via a short bridging group) to a sphingoid/ceramide backbone"


# Example usage (for testing purposes; uncomment to run):
# smiles_example = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_glycosphingolipid(smiles_example)
# print(result, reason)