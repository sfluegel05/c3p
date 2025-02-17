"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid

A glycosphingolipid is defined as a glycolipid that is a carbohydrate‐containing derivative 
of a sphingoid or ceramide, in which the carbohydrate residue is attached (typically via an O‐
glycosidic linkage) to the O-1 position of the sphingoid.
"""

from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    
    Our improved heuristic strategy is:
      1. Identify candidate sugar rings by using RDKit’s ring info. A sugar ring is
         heuristically approximated as a five- or six-membered ring that contains exactly one oxygen atom.
      2. Identify a sphingoid/ceramide backbone by detecting a ceramide-like motif.
         Here we use the pattern "C(=O)NC(CO)" (ignoring explicit chirality) which captures many ceramide fragments.
      3. Look for a glycosidic linkage between one of the sugar rings and the sphingoid fragment.
         Instead of checking only for a directly attached oxygen, we also allow a bridging atom (e.g. phosphorus).
      4. To avoid some false positives, we require that the molecule has a minimum number of carbons.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a glycosphingolipid, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Optional: Check that the molecule is large enough to be a lipid.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons to be a glycosphingolipid"
    
    # Step 1: Identify candidate sugar rings.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # each ring is recorded as a set of atom indices
    for ring in ring_info:
        if len(ring) in (5, 6):
            # For a typical pyranose or furanose, expect exactly one oxygen atom in the ring.
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count == 1:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety detected"
    
    # Step 2: Identify sphingoid/ceramide backbone.
    # Use a relaxed SMARTS for a ceramide-like pattern.
    sphingoid_smarts = "C(=O)NC(CO)"
    sphingo_pattern = Chem.MolFromSmarts(sphingoid_smarts)
    sphingo_matches = mol.GetSubstructMatches(sphingo_pattern)
    if not sphingo_matches:
        return False, "No sphingoid/ceramide backbone detected"
    # Collect all atom indices in any sphingoid match.
    sphingo_atoms = set()
    for match in sphingo_matches:
        sphingo_atoms.update(match)
    
    # Step 3: Look for a glycosidic linkage.
    # For each sugar ring candidate, we search for a ring carbon (candidate anomeric carbon)
    # that is attached to an exocyclic oxygen. We then check whether that oxygen is directly bonded
    # (or via a bridging phosphorus atom) to an atom from the sphingoid motif.
    glyco_link_found = False
    for ring in sugar_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Consider only carbon atoms in the sugar ring.
            if atom.GetAtomicNum() != 6:
                continue
            # Iterate over neighbors not in the ring.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                # Look for an exocyclic oxygen.
                if nbr.GetAtomicNum() == 8:
                    # Check that this oxygen has a connectivity suggesting a bridging group (often only two heavy atoms).
                    heavy_neighbors = [n for n in nbr.GetNeighbors() if n.GetAtomicNum() > 1]
                    if len(heavy_neighbors) != 2:
                        continue
                    # Option A: Direct linkage.
                    for n2 in nbr.GetNeighbors():
                        if n2.GetIdx() != atom_idx and n2.GetIdx() in sphingo_atoms:
                            glyco_link_found = True
                            break
                    if glyco_link_found:
                        break
                    # Option B: Allow for bridging via phosphorus.
                    for n2 in nbr.GetNeighbors():
                        if n2.GetAtomicNum() == 15:  # phosphorus
                            for n3 in n2.GetNeighbors():
                                if n3.GetIdx() not in (nbr.GetIdx(), atom_idx) and n3.GetIdx() in sphingo_atoms:
                                    glyco_link_found = True
                                    break
                        if glyco_link_found:
                            break
                if glyco_link_found:
                    break
            if glyco_link_found:
                break
        if glyco_link_found:
            break

    if not glyco_link_found:
        return False, "Sugar and sphingoid moieties are not glycosidically linked (even allowing one bridging atom)"
    
    return True, "Molecule contains a glycosphingolipid structure: a sugar ring is glycosidically linked (directly or via a bridging phosphorus) to a sphingoid/ceramide backbone"

# Example usage (for testing; uncomment to run):
# smiles_example = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_glycosphingolipid(smiles_example)
# print(result, reason)