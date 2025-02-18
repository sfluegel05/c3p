"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition:
  Any member of the class of chlorobenzenes carrying exactly four chlorine atoms attached directly 
  to a benzene ring. Additional substituents are tolerated provided they are “simple” – that is, 
  they are either a single atom substituent (like –OH, –CN, –F, etc.) or a biphenyl linkage 
  where the attached ring is itself a clean benzene ring.
  
Improvement over the previous version:
  • When looking at a candidate six-membered aromatic (benzene) ring, for every substituent attached 
    to a ring atom (neighbors not in the ring) we count an atom if it is chlorine; otherwise we 
    examine the connected substituent fragment.
  • The fragment is “simple” if (a) it has three or fewer heavy atoms or (b) it is exactly 6 heavy atoms 
    and all atoms are aromatic carbons (a benzene ring).
  • If any substituent is not “simple”, that candidate ring is disqualified.
  • If more than one ring is present we first sanitize the molecule and then focus on the largest fragment.
  
If the analysis is too ambiguous the function returns (None, None).
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule belongs to the tetrachlorobenzene class.
    
    It examines every six-membered aromatic ring (a benzene candidate) in the largest fragment 
    of the molecule. For each candidate ring it counts chlorine atoms directly attached to ring carbons.
    Other substituents are allowed only if they are “simple” – a small fragment (3 or fewer heavy atoms)
    or a direct benzene linkage (6 heavy atoms all aromatic carbons). If a candidate ring has exactly 4 
    directly bonded Cl atoms and every extra substituent is acceptable then the molecule is classified as a 
    tetrachlorobenzene. Otherwise it is not.
    
    Args:
      smiles (str): SMILES string of the molecule.
    
    Returns:
      (bool, str): Tuple where the first element is True if a qualifying benzene ring is found,
                   False if not, and (None, None) if the analysis was ambiguous.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Sanitize molecule
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {e}"
    
    # If multiple fragments exist, work on the largest one
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())

    ring_info = mol.GetRingInfo()
    if not ring_info or ring_info.NumRings() == 0:
        return False, "No rings found in the molecule"

    # Helper: perform a DFS starting from atom (neigh) but do not cross any atoms in blocked_set.
    def get_fragment_atoms(start_atom, blocked_set):
        to_visit = [start_atom.GetIdx()]
        visited = set()
        while to_visit:
            current = to_visit.pop()
            if current in visited:
                continue
            visited.add(current)
            atom = mol.GetAtomWithIdx(current)
            for nb in atom.GetNeighbors():
                if nb.GetIdx() in blocked_set:
                    continue
                if nb.GetIdx() not in visited:
                    to_visit.append(nb.GetIdx())
        return visited

    # Helper: determine if the substituent on the candidate ring (starting from neigh)
    # is "simple": either containing 3 or fewer heavy atoms, or exactly 6 heavy atoms that form a benzene ring.
    def is_allowed_substituent(neigh, candidate_idxs):
        # Reject if formal charge exists.
        if neigh.GetFormalCharge() != 0:
            return False
        # Use DFS to get the fragment (exclude candidate ring atoms)
        frag_atoms = get_fragment_atoms(neigh, candidate_idxs)
        # Count heavy atoms (atomic number > 1)
        heavy_atoms = [mol.GetAtomWithIdx(idx) for idx in frag_atoms if mol.GetAtomWithIdx(idx).GetAtomicNum() > 1]
        count = len(heavy_atoms)
        if count <= 3:
            return True
        if count == 6:
            # Check if the fragment is a benzene: every atom must be carbon and aromatic.
            if all(atom.GetAtomicNum() == 6 and atom.GetIsAromatic() for atom in heavy_atoms):
                return True
        return False

    # Loop over candidate rings: only examine rings that are 6-membered and all aromatic carbons.
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

        cl_count = 0  # number of directly attached Cl atoms
        disallowed_extras = 0  # count of substituents that are not allowed
        # Iterate over atoms in the candidate benzene ring.
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            for neigh in atom.GetNeighbors():
                # Only consider atoms not in the candidate ring.
                if neigh.GetIdx() in candidate_ring:
                    continue
                # If substituent is chlorine (atomic number 17), count as Cl.
                if neigh.GetAtomicNum() == 17:
                    cl_count += 1
                    continue
                # Otherwise, check if the attached substituent is allowed.
                if not is_allowed_substituent(neigh, set(candidate_ring)):
                    disallowed_extras += 1
        # Check if candidate benzene ring qualifies as tetrachlorobenzene:
        if cl_count == 4 and disallowed_extras == 0:
            return True, ("Found benzene ring with exactly 4 chlorine substituents and "
                          "all extra substituents are simple (small fragment or benzene linkage)")
    
    # If no candidate qualifies, return a failure reason.
    return False, "No benzene ring with exactly 4 chlorine substituents and acceptable extra groups found"


# Example usage:
if __name__ == "__main__":
    # List a few test SMILES for demonstration:
    test_cases = [
        ("Clc1cc(Cl)c(Cl)c(Cl)c1", "1,2,3,5-tetrachlorobenzene"),
        ("Clc1ccc(Cl)c(Cl)c1Cl", "1,2,3,4-tetrachlorobenzene"),
        ("Oc1c(Cl)c(Cl)cc(Cl)c1Cl", "2,3,5,6-tetrachlorophenol"),
        ("ClC=1C(C=2C(Cl)=CC(Cl)=C(Cl)C2)=CC(Cl)=C(Cl)C1Cl", "PCB180 (false positive example?)"),
    ]
    for smi, name in test_cases:
        result, reason = is_tetrachlorobenzene(smi)
        print(f"{name}: {result} -- {reason}")