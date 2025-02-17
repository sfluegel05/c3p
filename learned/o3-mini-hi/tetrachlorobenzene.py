"""
Classifies: CHEBI:26888 tetrachlorobenzene
"""
"""
Classifies: Tetrachlorobenzene
Definition:
  Any member of the class of chlorobenzenes carrying exactly four chlorine atoms attached directly 
  to a benzene ring. Additional substituents are tolerated provided they are “simple” – that is, 
  either a single atom substituent (like –OH, –CN, –F, etc.) or a biphenyl linkage where the attached 
  ring is itself a clean benzene ring.
"""

from rdkit import Chem

def is_tetrachlorobenzene(smiles: str):
    """
    Determines if a molecule belongs to the tetrachlorobenzene class.
    
    It considers every six-membered aromatic (benzene) ring in the largest fragment of the molecule.
    For each candidate ring, it counts chlorine atoms (atomic number 17) directly attached to the ring.
    For any other substituent, it performs a DFS (excluding ring atoms) to obtain the fragment 
    and checks if the substituent is "simple" – meaning 3 or fewer heavy atoms, or exactly 6 heavy atoms
    that form a benzene ring (all aromatic carbons). Only when the candidate ring has exactly 4 chlorine 
    substituents and no unacceptable extras is the molecule classified as tetrachlorobenzene.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple. The first element is True if a qualifying benzene ring is found, 
                   False if not, or (None, None) if the analysis is ambiguous.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Sanitize the molecule.
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization error: {e}"
    
    # If multiple fragments exist, focus on the largest fragment.
    frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    if len(frags) > 1:
        mol = max(frags, key=lambda m: m.GetNumHeavyAtoms())
    
    ring_info = mol.GetRingInfo()
    if not ring_info or ring_info.NumRings() == 0:
        return False, "No rings found in the molecule"
    
    # Helper function: perform a DFS starting from an atom (start_atom) but do not cross atoms
    # that are in the blocked_set (e.g., the candidate ring atoms)
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

    # Helper function: determine if the substituent fragment starting from neighbor "neigh"
    # (with candidate ring atoms blocked) is "simple":
    #   i) it has 3 or fewer heavy atoms (atomic num > 1)
    #   ii) OR it has exactly 6 heavy atoms and each is an aromatic carbon (benzene ring fragment)
    def is_allowed_substituent(neigh, candidate_idxs):
        # Reject if formal charge exists (could complicate definitions)
        if neigh.GetFormalCharge() != 0:
            return False
        frag_atom_idxs = get_fragment_atoms(neigh, candidate_idxs)
        # Determine heavy atoms (atomic number > 1)
        heavy_atoms = [mol.GetAtomWithIdx(i) for i in frag_atom_idxs if mol.GetAtomWithIdx(i).GetAtomicNum() > 1]
        count = len(heavy_atoms)
        if count <= 3:
            return True
        if count == 6:
            # Check if every heavy atom is carbon and aromatic (i.e. the fragment is benzene)
            if all(a.GetAtomicNum() == 6 and a.GetIsAromatic() for a in heavy_atoms):
                return True
        return False

    # Loop over candidate rings: consider only rings of 6 atoms that are all aromatic carbons.
    for candidate_ring in ring_info.AtomRings():
        if len(candidate_ring) != 6:
            continue
        is_candidate_benzene = True
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            if atom.GetAtomicNum() != 6 or not atom.GetIsAromatic():
                is_candidate_benzene = False
                break
        if not is_candidate_benzene:
            continue
        
        cl_count = 0  # count of chlorine substituents directly attached
        disallowed_extras = 0  # count of substituents which are not allowed (i.e. not "simple")
        # Iterate over atoms in the candidate benzene ring.
        for idx in candidate_ring:
            atom = mol.GetAtomWithIdx(idx)
            for neigh in atom.GetNeighbors():
                # Process substituents only if Neigh is not in the candidate ring.
                if neigh.GetIdx() in candidate_ring:
                    continue
                # If the substituent atom is chlorine, count it.
                if neigh.GetAtomicNum() == 17:
                    cl_count += 1
                else:
                    # For any other substituent, check if the fragment is allowed.
                    if not is_allowed_substituent(neigh, set(candidate_ring)):
                        disallowed_extras += 1
        # A candidate ring qualifies if it has exactly 4 chlorine substituents and no disallowed extras.
        if cl_count == 4 and disallowed_extras == 0:
            return True, ("Found benzene ring with exactly 4 chlorine substituents and "
                          "all extra substituents are simple (small fragment or benzene linkage)")
    
    return False, "No benzene ring with exactly 4 chlorine substituents and acceptable extra groups found"


# Example usage (only for testing - can be removed in production):
if __name__ == "__main__":
    # A set of test cases including both true positives and negatives based on provided outcomes:
    test_cases = [
        ("Clc1cc(Cl)c(Cl)c(Cl)c1", "1,2,3,5-tetrachlorobenzene (TP)"),
        ("Clc1ccc(Cl)c(Cl)c1Cl", "1,2,3,4-tetrachlorobenzene (TP)"),
        ("Oc1c(Cl)c(Cl)cc(Cl)c1Cl", "2,3,5,6-tetrachlorophenol (TP)"),
        ("C1(O)=C(C(=C(C=C1Cl)Cl)Cl)Cl", "2,3,4,5-tetrachlorophenol (TP)"),
        ("Clc1c(Cl)c(C#N)c(Cl)c(C#N)c1Cl", "chlorothalonil (TP)"),
        ("Oc1c(O)c(Cl)c(Cl)c(Cl)c1Cl", "tetrachlorocatechol (TP)"),
        ("Clc1cc(Cl)c(Cl)c(-c2ccccc2)c1Cl", "2,3,5,6-tetrachlorobiphenyl (TP)"),
        ("Clc1cc(-c2ccccc2)c(Cl)c(Cl)c1Cl", "2,3,4,5-tetrachlorobiphenyl (TP)"),
        ("Clc1cc(Cl)c(Cl)cc1Cl", "1,2,4,5-tetrachlorobenzene (TP)"),
        ("Oc1c(Cl)c(Cl)c(O)c(Cl)c1Cl", "tetrachlorohydroquinone (TP)"),
        ("C1(O)=C(C(=C(C=C1Cl)Cl)Cl)Cl", "2,3,4,6-tetrachlorophenol (TP)"),
        # False positives / negatives from previous attempt:
        ("ClC1=C(Cl)C(Cl)=C(O)C(=C1Cl)OC", "Tetrachloropyrocatechol methyl ether (FP)"),
        ("ClC1=C(OC)C(Cl)=CC(Cl)=C1Cl", "1,2,3,5-Tetrachloro-4-methoxybenzene (FP)"),
        ("ClC1=C(C(Cl)=C(Cl)C(=C1Cl)C(O)=O)C(O)=O", "Chlorthal (FP)"),
        ("Clc1cccc(c1)-c1cc(Cl)c(Cl)c(Cl)c1Cl", "2,3,3',4,5-pentachlorobiphenyl (FN)"),
    ]
    
    for smi, label in test_cases:
        result, reason = is_tetrachlorobenzene(smi)
        print(f"{label}: {result} -- {reason}")