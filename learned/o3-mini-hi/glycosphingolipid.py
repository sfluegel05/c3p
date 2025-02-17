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
         heuristically approximated as a five- or six-membered ring containing exactly one
         ring oxygen.
      2. Identify a sphingoid/ceramide backbone by detecting a ceramide-like motif.
         Here we use a relaxed SMARTS "C(=O)NC(CO)" and collect all atoms from any match.
      3. Look for a glycosidic linkage between one of the sugar rings and the sphingoid fragment.
         Instead of forcing the exocyclic oxygen to have exactly two heavy neighbors, we simply:
           • For each carbon in the sugar ring (a candidate anomeric carbon) that has at least one
             exocyclic oxygen substituent (i.e. a neighbor not in the ring that is oxygen),
           • Check if that oxygen is directly bonded to an atom from the sphingoid motif or via
             a bridging phosphorus atom.
      4. To reduce false positives, we also require that the molecule contains at least 20 carbon atoms.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is classified as a glycosphingolipid, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is large enough to be a lipid.
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 20:
        return False, "Too few carbons to be a glycosphingolipid"
    
    # Step 1: Identify candidate sugar rings.
    # Heuristically, a sugar ring is 5- or 6-membered and contains exactly one ring oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # list of sets of atom indices
    for ring in ring_info:
        if len(ring) in (5, 6):
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count == 1:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety detected"
    
    # Step 2: Identify the sphingoid/ceramide backbone.
    # This SMARTS is a relaxed pattern for a ceramide-like fragment.
    sphingoid_smarts = "C(=O)NC(CO)"
    sphingo_pattern = Chem.MolFromSmarts(sphingoid_smarts)
    sphingo_matches = mol.GetSubstructMatches(sphingo_pattern)
    if not sphingo_matches:
        return False, "No sphingoid/ceramide backbone detected"
    # Collect all atom indices from any sphingoid match.
    sphingo_atoms = set()
    for match in sphingo_matches:
        sphingo_atoms.update(match)
    
    # Step 3: Look for a glycosidic linkage.
    # For each sugar ring candidate, look for a candidate anomeric carbon.
    # We assume an anomeric carbon is a carbon in the ring that has at least one exocyclic oxygen.
    glyco_link_found = False
    for ring in sugar_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Consider only carbon atoms.
            if atom.GetAtomicNum() != 6:
                continue
            # Find exocyclic oxygen substituents.
            for nbr in atom.GetNeighbors():
                if nbr.GetIdx() in ring:
                    continue
                if nbr.GetAtomicNum() == 8:
                    # We have a candidate exocyclic oxygen attached to this carbon.
                    # Option A: Direct linkage between this oxygen and the sphingoid backbone.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() != atom_idx and nbr2.GetIdx() in sphingo_atoms:
                            glyco_link_found = True
                            break
                    if glyco_link_found:
                        break
                    # Option B: Allow a bridging phosphorus.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() == atom_idx:
                            continue
                        if nbr2.GetAtomicNum() == 15:  # phosphorus
                            for nbr3 in nbr2.GetNeighbors():
                                if nbr3.GetIdx() not in (nbr.GetIdx(), atom_idx) and nbr3.GetIdx() in sphingo_atoms:
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

# Example usage (for testing purposes; uncomment to run):
# smiles_example = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_glycosphingolipid(smiles_example)
# print(result, reason)