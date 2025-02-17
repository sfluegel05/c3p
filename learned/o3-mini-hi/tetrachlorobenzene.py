"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition:
  Any member of the class of chlorobenzenes carrying four chlorine atoms attached directly 
  to a benzene ring. Valid molecules may have an extra substituent if it is either a simple group 
  (e.g. –OH or –CN) or a benzene ring linkage as in biphenyls. Extra substituents that carry exotic 
  features (e.g. nonzero formal charge, or larger multi–atom groups that are not directly a benzene ring)
  disfavor the classification.
  
Improved heuristic:
  1. Parse the molecule and (if needed) select its largest fragment.
  2. For every 6–membered ring composed of aromatic carbons (a benzene candidate) do:
      – For each ring atom, consider every neighbor that is not part of the ring.
      – If the neighbor is Cl, count it.
      – Otherwise, if the substituent is allowed (see allowed_substituent) then ignore it.
         Otherwise, flag the ring as having disallowed extra groups.
      – If the candidate benzene ring has exactly four Cl substituents and no disallowed extras, 
        then return True (with a proper explanation).
  3. If no candidate ring qualifies, return False.
  
Note: This heuristic is not perfect – if the analysis is too ambiguous the function returns (None, None).
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule is a tetrachlorobenzene based on its SMILES string.
  
    The function first looks for candidate benzene rings (6 aromatic carbons). For each candidate, it
    counts how many chlorine atoms are directly bonded to ring carbons. Non–chlorine substituents
    are allowed only if they are simple groups – either a single atom such as –OH, a –CN group,
    or a direct benzene–benzene (biphenyl) linkage.
  
    Args:
       smiles (str): SMILES string of the molecule.
  
    Returns:
       (bool, str): A tuple where the first element is True if a qualifying benzene ring is found,
                    and False otherwise; the second element gives a reason.
       If the analysis cannot be confidently done, returns (None, None).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {e}"

    # If multiple fragments are present, use the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    ring_info = mol.GetRingInfo()
    if not ring_info or ring_info.NumRings() == 0:
        return False, "No rings found in the molecule"
    
    # Helper: decide whether a substituent (atom) is allowed.
    # Allowed if:
    #   - It is chlorine (handled separately in the main loop).
    #   - It is a small simple group such as –OH (oxygen with degree 1).
    #   - It is a carbon that is part of a cyanide (triple bond to nitrogen).
    #   - It is an aromatic carbon that belongs to a six-membered aromatic ring
    #     (which we interpret as a biphenyl linkage).
    def allowed_substituent(neigh, candidate_ring):
        # Exclude substituents carrying a formal charge.
        if neigh.GetFormalCharge() != 0:
            return False
        symbol = neigh.GetSymbol()
        # Allow single–atom substituents (–OH, –F, –I, etc.) except we already count Cl separately.
        if neigh.GetDegree() == 1:
            return True
        # Allow a carbon that is connected via a triple bond to a nitrogen (–CN)
        if symbol == "C":
            for bond in neigh.GetBonds():
                if bond.GetBondType() == Chem.BondType.TRIPLE:
                    other = bond.GetOtherAtom(neigh)
                    if other.GetAtomicNum() == 7:
                        return True
        # Allow biphenyl linkages: if the neighbor is an aromatic carbon that belongs to a six-membered ring
        # (other than the candidate ring), then permit.
        if symbol == "C" and neigh.GetIsAromatic():
            for r in ring_info.AtomRings():
                if len(r) == 6 and neigh.GetIdx() in r:
                    # check that this ring is not the candidate ring
                    if set(r) != set(candidate_ring):
                        return True
        return False

    # Loop over rings: only consider 6-membered rings that are all aromatic carbons (benzene candidates).
    for candidate_ring in ring_info.AtomRings():
        if len(candidate_ring) != 6:
            continue
        is_benzene = True
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_benzene = False
                break
        if not is_benzene:
            continue
        
        # Examine substituents attached directly to the benzene candidate.
        cl_count = 0
        disallowed_extras = 0
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            for neigh in atom.GetNeighbors():
                # Only consider atoms that are not in the candidate ring.
                if neigh.GetIdx() in candidate_ring:
                    continue
                if neigh.GetAtomicNum() == 17:
                    cl_count += 1
                else:
                    if not allowed_substituent(neigh, candidate_ring):
                        disallowed_extras += 1
        # To qualify, require exactly 4 Cl atoms and no disallowed extra substituents.
        if cl_count == 4 and disallowed_extras == 0:
            return True, "Found benzene ring with exactly 4 chlorine substituents and acceptable additional groups"
    
    # If none of the candidate benzene rings qualifies, report failure.
    return False, "No benzene ring with exactly 4 chlorine substituents (and acceptable extra groups) found"

# Example usage:
if __name__ == "__main__":
    # Test with one known tetrachlorobenzene:
    test_smiles = "Clc1cc(Cl)c(Cl)c(Cl)c1"
    is_tc, reason = is_tetrachlorobenzene(test_smiles)
    print(is_tc, reason)