"""
Classifies: CHEBI:48039 dihydroflavonols
"""
#!/usr/bin/env python3
"""
Classifies: dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.
This implementation works by locating a carbonyl (C=O) that is directly attached to a saturated carbon bearing an -OH.
Then it checks that these two atoms are part of a six-membered ring that also contains a heterocyclic oxygen.
Finally it verifies that a different carbon in that same ring (corresponding to C2 of flavanones)
bears an external aromatic substituent (the B‐ring).
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a given SMILES string corresponds to a dihydroflavonol.
    The key flavanone (2-phenylchroman-4-one) scaffold of dihydroflavonols is identified by:
      - a carbonyl group (C=O) attached to a saturated carbon that bears an –OH (expected to be C3),
      - both atoms are in a six-membered ring which also contains a ring oxygen (the heterocycle oxygen, expected to be at C1),
      - additionally, a different ring carbon (expected to be C2) is connected to an external aromatic group (the B‐ring).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule can be classified as a dihydroflavonol, False otherwise.
        str: Reason for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS to find a fragment: a carbonyl (C=O) directly bonded to a saturated carbon carrying an OH.
    # This is intended to capture the C4 (carbonyl) - C3 (with OH) fragment.
    query = Chem.MolFromSmarts("C(=O)[C]([OX2H])")
    if query is None:
        return None, None  # Should not happen.
    
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Molecule does not contain the key fragment C(=O)[C]([OX2H]) expected for dihydroflavonols"
    
    # Get all rings as sets of atom indices in the molecule.
    rings = [set(r) for r in Chem.GetSymmSSSR(mol)]
    
    # Loop over all fragment matches. For each, we expect:
    #   match[0] = carbonyl carbon (C4)
    #   match[1] = saturated carbon bearing -OH (C3)
    for match in matches:
        carbonyl_idx = match[0]
        c3_idx = match[1]
        
        # Look for a ring that is six members long and includes both carbonyl and c3.
        candidate_rings = [ring for ring in rings if carbonyl_idx in ring and c3_idx in ring and len(ring) == 6]
        if not candidate_rings:
            continue  # Try next match
        
        # For each such ring, check for two features:
        #   (i) The ring must contain a heterocyclic oxygen (expected to be the ring oxygen at C1).
        #   (ii) There should be a ring carbon (expected C2) that is connected to an external aromatic group.
        for ring in candidate_rings:
            # (i) Find a ring oxygen atom.
            ring_oxygens = [idx for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 
                             and mol.GetAtomWithIdx(idx).GetSymbol() == "O"]
            if not ring_oxygens:
                continue  # No ring oxygen; skip this ring
            
            # (ii) Look for a candidate ring carbon that is attached to an aromatic substituent outside the ring.
            aromatic_attachment_found = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                # Consider only carbon atoms that are part of the ring and are not the carbonyl atom or the C3 (with -OH).
                if atom.GetAtomicNum() != 6 or idx in (carbonyl_idx, c3_idx):
                    continue
                # To be a candidate for C2, it should be bonded to the ring oxygen.
                neighbors_in_ring = set(nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetIdx() in ring)
                if not neighbors_in_ring.intersection(ring_oxygens):
                    continue  # Not adjacent to ring oxygen
                # Now, check if this atom has at least one aromatic neighbor that is not part of the same ring.
                for nbr in atom.GetNeighbors():
                    if nbr.GetIdx() not in ring and nbr.GetIsAromatic():
                        aromatic_attachment_found = True
                        break
                if aromatic_attachment_found:
                    break  # We found our candidate (expected C2 with the B‐ring)
            
            # If both criteria (ring oxygen and aromatic attachment at a different ring carbon) are satisfied, classify as dihydroflavonol.
            if aromatic_attachment_found:
                return True, ("Matches dihydroflavonol core: contains a 2-phenylchroman-4-one (flavanone) scaffold with a -OH at the C3 position.")
        
    return False, "Molecule does not match the dihydroflavonol core criteria."

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_examples = [
        # Provided examples
        "O1C(C(O)C(=O)C=2C1=CC=3OCOC3C2O",  # 3,5-Dihydroxy-6,7-methylenedioxyflavanone
        "O[C@@H]1[C@H](Oc2cc(O)ccc2C1=O)c1ccc(O)cc1",    # garbanzol
        "CC(C)=CCc1c(O)cc(O)c2C(=O)[C@H](O)[C@H](Oc12)c1ccccc1",  # glepidotin B
        "O1[C@H]([C@H](O)C(=O)C=2C1=CC(O)=CC2O)C3=CC=C(O)C=C3",  # (2S,3S)-3,5,7-Trihydroxy-2-(4-hydroxyphenyl)-2,3-dihydro-4H-chromen-4-one
        "COC1=C(O)C=CC(=C1)[C@H]1OC2=CC(O)=CC(O)=C2C(=O)[C@@H]1O", # dihydroisorhamnetin
        "OC1C(Oc2cc(O)cc(O)c2C1=O)c1ccc(O)c(O)c1",           # (+)-taxifolin
        "O[C@H]1[C@@H](Oc2cc(O)cc(O)c2C1=O)c1cc(O)c(O)c(O)c1",  # (-)-dihydromyricetin
    ]
    for s in test_examples:
        result, reason = is_dihydroflavonols(s)
        print(f"SMILES: {s}\nClassified as dihydroflavonol? {result}\nReason: {reason}\n")