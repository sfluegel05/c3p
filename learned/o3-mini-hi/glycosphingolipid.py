"""
Classifies: CHEBI:24402 glycosphingolipid
"""
"""
Classifies: Glycosphingolipid
A glycosphingolipid is defined as a glycolipid that is a carbohydrate‐containing derivative 
of a sphingoid or ceramide, with the carbohydrate residue attached by a glycosidic linkage 
to O-1 of the sphingoid.
"""

from rdkit import Chem

def is_glycosphingolipid(smiles: str):
    """
    Determines if a molecule is a glycosphingolipid based on its SMILES string.
    
    Our heuristic strategy is:
      1. Identify a carbohydrate ring by iterating over all rings and flagging those 
         that are five- or six-membered and contain exactly one oxygen atom.
      2. Identify a sphingoid/ceramide backbone by detecting a fragment 
         "C(=O)NC(CO)" which heuristically captures many ceramide-like motifs.
      3. Verify that one of the sugar rings is connected via an exocyclic oxygen 
         (i.e. a glycosidic linkage) to the sphingoid fragment.
         
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if the molecule is a glycosphingolipid, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 1: Detect candidate sugar rings.
    # We use ring info from RDKit and select rings of 5 or 6 atoms that contain exactly one oxygen.
    ring_info = mol.GetRingInfo().AtomRings()
    sugar_rings = []  # list of sets; each set contains the indices of atoms in a candidate sugar ring
    for ring in ring_info:
        if len(ring) in (5, 6):
            # Count oxygen atoms inside the ring
            oxy_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8)
            if oxy_count == 1:
                sugar_rings.append(set(ring))
    if not sugar_rings:
        return False, "No carbohydrate (sugar) moiety detected"

    # Step 2: Detect a sphingoid/ceramide backbone.
    # Here we use a pattern that captures an amide group attached to a carbon that carries a -CH2OH.
    sphingoid_smarts = "C(=O)NC(CO)"
    sphingo_pattern = Chem.MolFromSmarts(sphingoid_smarts)
    sphingo_matches = mol.GetSubstructMatches(sphingo_pattern)
    if not sphingo_matches:
        return False, "No sphingoid/ceramide backbone detected"
    # Collect all atom indices involved in any sphingoid match
    sphingo_atoms = set()
    for match in sphingo_matches:
        sphingo_atoms.update(match)
    
    # Step 3: Check for a glycosidic linkage.
    # For each candidate sugar ring, check whether one of its carbon atoms is bonded to an exocyclic oxygen
    # that in turn is bonded to at least one atom from the sphingoid fragment.
    glyco_link_found = False
    for ring in sugar_rings:
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Typically, the anomeric carbon is a carbon in the sugar; check only if this atom is carbon.
            if atom.GetAtomicNum() != 6:
                continue
            # Look at neighbors not in this ring; such oxygen(s) (if any) could participate in a glycosidic bond.
            for nbr in atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8 and nbr.GetIdx() not in ring:
                    # Check if this oxygen also bonds to an atom in the sphingoid fragment.
                    for nbr2 in nbr.GetNeighbors():
                        if nbr2.GetIdx() != atom_idx and nbr2.GetIdx() in sphingo_atoms:
                            glyco_link_found = True
                            break
                    if glyco_link_found:
                        break
            if glyco_link_found:
                break
        if glyco_link_found:
            break

    if not glyco_link_found:
        return False, "Sugar and sphingoid moieties are not glycosidically linked"

    return True, "Molecule contains a glycosphingolipid structure: carbohydrate moiety glycosidically linked to a sphingoid/ceramide backbone"

# Example usage (for testing; uncomment to run):
# smiles_example = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO[C@@H]1O[C@H](CO)[C@@H](O)[C@H](O)[C@H]1O)[C@H](O)\\C=C\\CCCCCCCCCCCCC"
# result, reason = is_glycosphingolipid(smiles_example)
# print(result, reason)