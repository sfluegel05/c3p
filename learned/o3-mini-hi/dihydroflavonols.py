"""
Classifies: CHEBI:48039 dihydroflavonols
"""
"""
Classifies: dihydroflavonols
Definition: Any hydroxyflavanone in which a hydroxy group is present at position 3 of the heterocyclic ring.
This implementation does not rely on a single rigid SMARTS but uses a fragment search combined with ring analysis.
Key criteria applied:
  - The molecule must contain a fragment “C(=O)[CX4]([OX2H])”,
    meaning that a carbonyl (C=O) is bonded to a saturated carbon bearing an –OH group.
  - These two atoms must belong to a six-membered ring.
  - The ring must contain a heterocyclic oxygen.
  - The saturated carbon (bearing –OH) should also be bonded (outside the ring) to an aromatic group (the B‐ring).
"""

from rdkit import Chem

def is_dihydroflavonols(smiles: str):
    """
    Determines if a given SMILES string corresponds to a dihydroflavonol.
    
    This function first searches for the fragment "C(=O)[CX4]([OX2H])" that
    ensures a carbonyl group is directly bonded to a saturated carbon with a hydroxyl.
    For each match found, it verifies that:
      - Both atoms are part of the same ring that has 6 members.
      - That ring contains at least one oxygen atom (ensuring a fused benzopyran).
      - The carbon with the –OH has at least one neighbor (outside the ring)
        that is aromatic (expected to be the B‐ring).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a dihydroflavonol, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for a carbonyl (C=O) attached to a tetrahedral carbon bearing an -OH.
    query = Chem.MolFromSmarts("C(=O)[CX4]([OX2H])")
    if query is None:
        return None, None  # Should not happen.
    
    matches = mol.GetSubstructMatches(query)
    if not matches:
        return False, "Molecule does not contain the key fragment C(=O)[CX4]([OX2H]) required for dihydroflavonols"
    
    # Retrieve all rings (as tuples of atom indices) in the molecule.
    rings = [set(ring) for ring in Chem.GetSymmSSSR(mol)]
    
    # Loop over all matches found.
    for match in matches:
        # match[0] is the carbonyl carbon, match[1] is the saturated carbon with the OH (candidate C3).
        carbonyl_idx = match[0]
        c3_idx = match[1]
        
        # Check if these two atoms are in a six-membered ring that also contains an oxygen.
        found_ring = False
        for ring in rings:
            if carbonyl_idx in ring and c3_idx in ring and len(ring) == 6:
                # Check if at least one atom in the ring is oxygen.
                if any(mol.GetAtomWithIdx(idx).GetAtomicNum() == 8 and mol.GetAtomWithIdx(idx).GetSymbol() == "O" for idx in ring):
                    found_ring = True
                    break
        if not found_ring:
            continue  # Check next match
        
        # Check that the candidate C3 atom is connected to an aromatic atom not in the ring.
        c3_atom = mol.GetAtomWithIdx(c3_idx)
        aromatic_neighbor = False
        for nbr in c3_atom.GetNeighbors():
            if nbr.GetIsAromatic() and (nbr.GetIdx() not in ring):
                aromatic_neighbor = True
                break
        if not aromatic_neighbor:
            # Sometimes the aromatic B-ring might be fused into the system.
            # We check if there is any aromatic neighbor at all.
            for nbr in c3_atom.GetNeighbors():
                if nbr.GetIsAromatic():
                    aromatic_neighbor = True
                    break
                
        if not aromatic_neighbor:
            continue  # Not matching the expected B-ring connectivity
        
        # If all conditions are met, classify as dihydroflavonol.
        return True, ("Matches dihydroflavonol core: contains a benzopyran-4-one (flavanone) scaffold with a "
                      "saturated C2–C3 bond, a hydroxyl group at C3, and connectivity to an aromatic B‐ring.")
    
    return False, "Molecule does not match the dihydroflavonol core criteria."

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = [
        "O1C(C(O)C(=O)C=2C1=CC=3OCOC3C2O)C4=CC=CC=C4",  # 3,5-Dihydroxy-6,7-methylenedioxyflavanone
        "O[C@@H]1[C@H](Oc2cc(O)ccc2C1=O)c1ccc(O)cc1",    # garbanzol
        "CC(C)=CCc1c(O)cc(O)c2C(=O)[C@H](O)[C@H](Oc12)c1ccccc1",  # glepidotin B
        "O1[C@H]([C@H](O)C(=O)C=2C1=CC(O)=CC2O)C3=CC=C(O)C=C3",  # (2S,3S)-3,5,7-Trihydroxy-2-(4-hydroxyphenyl)-2,3-dihydro-4H-chromen-4-one
        "COC1=C(O)C=CC(=C1)[C@H]1OC2=CC(O)=CC(O)=C2C(=O)[C@@H]1O", # dihydroisorhamnetin
    ]
    for s in test_smiles:
        result, reason = is_dihydroflavonols(s)
        print(f"SMILES: {s}\nClassified as dihydroflavonol? {result}\nReason: {reason}\n")