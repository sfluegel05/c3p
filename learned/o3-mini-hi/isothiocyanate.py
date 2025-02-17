"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate (R-N=C=S)

An isothiocyanate is defined as an organosulfur compound with the general formula R-N=C=S.
This function inspects all occurrences of the N=C=S moiety and accepts the compound if at least one shows
a terminal, “free‐standing” N=C=S fragment (i.e. the -N=C=S group can be disconnected cleanly by breaking
the bond between the nitrogen and its substituent R).
"""

from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    For our purposes the molecule must contain at least one terminal group
    of the form R-N=C=S (where the -N=C=S unit is terminal, i.e. not embedded in a network).
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an isothiocyanate, False otherwise.
        str: Reason for the classification.
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # List of SMARTS patterns to identify possible isothiocyanate groups.
    patterns = ["N=C=S", "S=C=N"]
    candidate_matches = []
    
    # Get all substructure matches for each pattern.
    for pat in patterns:
        pat_mol = Chem.MolFromSmarts(pat)
        if not pat_mol:
            continue
        matches = mol.GetSubstructMatches(pat_mol)
        candidate_matches.extend(matches)
    
    # Deduplicate matches by sorting the tuple of indices.
    unique_matches = []
    seen = set()
    for match in candidate_matches:
        key = tuple(sorted(match))
        if key not in seen:
            seen.add(key)
            unique_matches.append(match)
    
    if not unique_matches:
        return False, "No isothiocyanate SMARTS match found"
    
    # Iterate over each unique candidate match.
    # We allow a molecule to qualify if at least one occurrence of a terminal isothiocyanate group is found.
    for match in unique_matches:
        # We expect 3 atoms. They might be in order (N, C, S) or (S, C, N). We'll reorder if needed.
        if len(match) != 3:
            continue  # unexpected match, skip
        first_atom = mol.GetAtomWithIdx(match[0])
        # If the first atom is sulfur (atomic num 16), assume the match came from "S=C=N" so reorder.
        if first_atom.GetAtomicNum() == 16:
            n_idx, c_idx, s_idx = match[2], match[1], match[0]
        else:
            n_idx, c_idx, s_idx = match[0], match[1], match[2]
        
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        s_atom = mol.GetAtomWithIdx(s_idx)
        
        # Check that bonds N-C and C-S exist and are double.
        bond_nc = mol.GetBondBetweenAtoms(n_idx, c_idx)
        bond_cs = mol.GetBondBetweenAtoms(c_idx, s_idx)
        if bond_nc is None or bond_cs is None:
            continue
        if bond_nc.GetBondType() != Chem.rdchem.BondType.DOUBLE or bond_cs.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            continue
        
        # Terminal connectivity check for the isothiocyanate group:
        # - Central carbon should have exactly 2 heavy neighbors (N and S).
        # - Sulfur must have exactly one heavy neighbor (the carbon).
        # - Nitrogen should have exactly 2 heavy neighbors; one is the carbon in the group and one is the substituent R.
        heavy_neighbors_c = [nbr.GetIdx() for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        heavy_neighbors_s = [nbr.GetIdx() for nbr in s_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        heavy_neighbors_n = [nbr.GetIdx() for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        
        if len(heavy_neighbors_c) != 2:
            continue
        if len(heavy_neighbors_s) != 1:
            continue
        if len(heavy_neighbors_n) != 2:
            continue
        
        # Identify the substituent R group: the neighbor of nitrogen NOT in the N=C=S group.
        r_idx = None
        for idx in heavy_neighbors_n:
            if idx != c_idx:
                r_idx = idx
                break
        if r_idx is None:
            continue
        
        # To ensure that the N=C=S group is terminal, we break the bond between the nitrogen and its R substituent
        # and then check if the group (N, C, S) becomes an isolated fragment.
        bond_nr = mol.GetBondBetweenAtoms(n_idx, r_idx)
        if bond_nr is None:
            continue
        bond_idx = bond_nr.GetIdx()
        
        # Fragment the molecule by breaking the N-R bond; add dummy atoms to mark the break.
        frag_mol = Chem.FragmentOnBonds(mol, [bond_idx], addDummies=True)
        frags = Chem.GetMolFrags(frag_mol, asMols=True, sanitizeFrags=True)
        
        # Now check if one of the fragments contains exactly the 3 atoms of our candidate.
        # Note that after fragmentation, dummy atoms (atomic num 0) may be present. We look for fragments whose heavy
        # atom indices (ignoring dummies) correspond exactly to {n_idx, c_idx, s_idx} from the original molecule.
        valid_terminal = False
        # First, we create a set of our candidate atom indices.
        candidate_set = {n_idx, c_idx, s_idx}
        for frag in frags:
            # Get the map numbers of atoms (if set during fragmentation) or try to retrieve original indices.
            atom_ids = []
            for atom in frag.GetAtoms():
                # Dummy atoms have atomicNum 0.
                if atom.GetAtomicNum() == 0:
                    continue
                # Use the atom property '_fragAtomMapping' if available;
                # if not, we can try to use GetIdx() but that index is local to the fragment.
                amap = atom.GetUnsignedProp("molAtomMapNumber") if atom.HasProp("molAtomMapNumber") else None
                if amap is not None:
                    atom_ids.append(amap)
            # Sometimes the fragmentation does not carry mapping numbers – in that case we fall back to matching size.
            # We compare fragment heavy atom count.
            if len(atom_ids) == 0:
                # if no mapping info is available, use number of heavy atoms.
                if frag.GetNumHeavyAtoms() == 3:
                    valid_terminal = True
                    break
            else:
                if set(atom_ids) == candidate_set:
                    valid_terminal = True
                    break
        
        if valid_terminal:
            return True, "Contains terminal isothiocyanate group (R-N=C=S)"
        
    # If no candidate match passes all tests, return False.
    return False, "No terminal isothiocyanate group (R-N=C=S) found"

# Example usage and testing
if __name__ == "__main__":
    # Test examples that should be classified as isothiocyanates.
    test_examples = {
        "sulforaphane": "N(=C=S)CCCCS(=O)C",
        "eosin 5-isothiocyanate": "Oc1c(Br)cc2c(Oc3c(Br)c(O)c(Br)cc3C22OC(=O)c3cc(ccc23)N=C=S)c1Br",
        "(2R)-2-hydroxy-3-butenyl isothiocyanate": "N(C[C@H](C=C)O)=C=S",
        "methyl isothiocyanate": "CN=C=S",
        "8-(methylsulfinyl)octyl isothiocyanate": "C(CCCCS(C)=O)CCCN=C=S",
        "isothiocyanato-Cyclohexane": "S=C=NC1CCCCC1",
        "(4-methoxyindol-3-ylmethyl)isothiocyanate": "COC1=CC=CC2=C1C(=CN2)CN=C=S",
        "1-Isothiocyanato-7-(methylthio)heptane": "S(CCCCCCCN=C=S)C",
        "3-hydroxypropyl isothiocyanate": "OCCCN=C=S",
        "Erucin": "S(CCCCN=C=S)C",
        "phenyl isothiocyanate": "S=C=Nc1ccccc1",
        "5-isothiocyanato-1-pentene": "C=CCCCN=C=S",
        "2-isothiocyanatoethyl phosphate": "S=C=NCCOP(O)(=O)O",
        "Isothiocyanatocyclopropane": "C1CC1N=C=S",
        "1-naphthyl isothiocyanate": "S=C=Nc1cccc2ccccc12",
        "(R)-sulforaphane": "N(=C=S)CCCC[S@](=O)C",
        "ethyl isothiocyanate": "CCN=C=S",
        "benzyl isothiocyanate": "S=C=NCc1ccccc1",
        "phenethyl isothiocyanate": "S=C=NCCc1ccccc1",
        "amoscanate": "C1=CC(=CC=C1NC=2C=CC(N=C=S)=CC2)[N+](=O)[O-]",
        "1-isothiocyanato-6-(methylsulfanyl)hexane": "S(CCCCCCN=C=S)C",
        "allyl isothiocyanate": "C=CCN=C=S",
        "toluene 2,6-diisothiocyanate": "Cc1c(cccc1N=C=S)N=C=S",
        "8-(methylthio)octylisothiocyanate": "S(CCCCCCCCN=C=S)C",
        "rhodamine B 5-isothiocyanate": "[Cl-].CCN(CC)c1ccc2c(-c3ccc(cc3C(O)=O)N=C=S)c3ccc(cc3oc2c1)=[N+](CC)CC",
        "Berteroin": "S(CCCCCN=C=S)C",
        "iberin": "N(CCCS(C)=O)=C=S",
        "Bis(4-isothiocyanatobutyl) disulfide": "S(SCCCCN=C=S)CCCCN=C=S",
        "toluene 2,4-diisothiocyanate": "Cc1ccc(cc1N=C=S)N=C=S",
    }
    
    for name, smi in test_examples.items():
        result, reason = is_isothiocyanate(smi)
        print(f"{name}: {result}, Reason: {reason}")