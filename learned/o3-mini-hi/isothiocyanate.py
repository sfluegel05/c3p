"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate (R-N=C=S)

An isothiocyanate is defined as an organosulfur compound with the general formula R-N=C=S.
This function looks for the SMARTS pattern in both directions, deduplicates duplicate matches,
and then enforces connectivity and exclusion checks.
"""

from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule is an isothiocyanate based on its SMILES string.
    Isothiocyanates are defined as having a terminal group R-N=C=S.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as an isothiocyanate, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for isothiocyanate in forward and reverse orders.
    patterns = ["N=C=S", "S=C=N"]
    candidate_matches = []
    for pat in patterns:
        pat_mol = Chem.MolFromSmarts(pat)
        if pat_mol is None:
            continue
        matches = mol.GetSubstructMatches(pat_mol)
        candidate_matches.extend(matches)
    
    # Deduplicate candidate matches by comparing sorted atom indices.
    unique_matches = []
    seen = set()
    for m in candidate_matches:
        key = tuple(sorted(m))
        if key not in seen:
            seen.add(key)
            unique_matches.append(m)
            
    # Reject if no candidate found.
    if not unique_matches:
        return False, "No isothiocyanate SMARTS match found"
    
    # Reject if more than one unique isothiocyanate moiety is found.
    if len(unique_matches) > 1:
        return False, f"Multiple isothiocyanate groups found ({len(unique_matches)} matches)"
    
    # Process the single candidate match.
    match = unique_matches[0]
    # By default assume match from pattern "N=C=S" gives order (N, C, S).
    n_idx, c_idx, s_idx = match[0], match[1], match[2]
    atom0 = mol.GetAtomWithIdx(match[0])
    # If the first atom is S then we likely got a match from "S=C=N" and we reorder.
    if atom0.GetAtomicNum() == 16:
        n_idx, c_idx, s_idx = match[2], match[1], match[0]
    
    n_atom = mol.GetAtomWithIdx(n_idx)
    c_atom = mol.GetAtomWithIdx(c_idx)
    s_atom = mol.GetAtomWithIdx(s_idx)
    
    # Check that the bonds N-C and C-S exist and are double bonds.
    bond_nc = mol.GetBondBetweenAtoms(n_idx, c_idx)
    bond_cs = mol.GetBondBetweenAtoms(c_idx, s_idx)
    if bond_nc is None or bond_cs is None:
        return False, "Expected bonds not found"
    if bond_nc.GetBondType() != Chem.rdchem.BondType.DOUBLE or bond_cs.GetBondType() != Chem.rdchem.BondType.DOUBLE:
        return False, "One or both bonds in the N=C=S group are not double bonds"
    
    # Enforce terminal heavy-atom connectivity.
    # Central carbon should have only N and S as heavy neighbors.
    heavy_neighbors_c = [nbr.GetIdx() for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if len(heavy_neighbors_c) != 2:
        return False, "Central carbon does not have exactly two heavy neighbors"
    
    # Terminal sulfur should have only the central carbon as heavy neighbor.
    heavy_neighbors_s = [nbr.GetIdx() for nbr in s_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if len(heavy_neighbors_s) != 1:
        return False, "Terminal sulfur is not terminal (unexpected heavy neighbors)"
    
    # Nitrogen should have exactly two heavy neighbors: one is the central carbon, and one is the R substituent.
    heavy_neighbors_n = [nbr.GetIdx() for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if len(heavy_neighbors_n) != 2:
        return False, "Nitrogen in the isothiocyanate group does not have exactly two heavy neighbors"
    
    # Identify the substituent R attached to nitrogen (the neighbor that is not the central C).
    r_idx = None
    for idx in heavy_neighbors_n:
        if idx != c_idx:
            r_idx = idx
            break
    if r_idx is None:
        return False, "Could not identify substituent (R) attached to nitrogen"
    r_atom = mol.GetAtomWithIdx(r_idx)
    
    # --- Exclusion filters on the R substituent ---
    # (1) If the R substituent is aromatic and directly attached to a carbonyl group, reject.
    if r_atom.GetIsAromatic():
        for nbr in r_atom.GetNeighbors():
            if nbr.GetIdx() == n_idx: 
                continue
            if nbr.GetAtomicNum() == 6:  # possible carbon that could be part of a carbonyl
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return False, "R substituent is aromatic and attached to a carbonyl group"
    
    # (2) If the R substituent (or its immediate neighbor) is linked to an oxidized heteroatom,
    # for instance S=O or P(=O), then reject.
    def has_oxidized_neighbor(atom):
        # Look for neighbor heteroatom with at least one double bond to oxygen.
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == n_idx:
                continue
            if nbr.GetAtomicNum() in (16, 15):  # S or P
                for bond in nbr.GetBonds():
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        other = bond.GetOtherAtom(nbr)
                        if other.GetAtomicNum() == 8:
                            return True
        return False

    if has_oxidized_neighbor(r_atom):
        return False, "R substituent is adjacent to an oxidized heteroatom (e.g. S=O or P=O)"
    for nbr in r_atom.GetNeighbors():
        if nbr.GetIdx() == n_idx:
            continue
        if has_oxidized_neighbor(nbr):
            return False, "Substituent branch near the isothiocyanate group is oxidized (e.g., S=O, P=O)"
    
    # (3) Enforce expected degrees of atoms in the group.
    # Ideal terminal group expected degrees: N:2, C:2, S:1.
    if n_atom.GetDegree() != 2 or c_atom.GetDegree() != 2 or s_atom.GetDegree() != 1:
        return False, "Degrees of atoms in N=C=S group are not as expected for a terminal isothiocyanate"
    
    return True, "Contains terminal isothiocyanate group (R-N=C=S)"

# Example usage and testing
if __name__ == "__main__":
    # Some test examples from the problem description.
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
    }
    
    for name, smi in test_examples.items():
        result, reason = is_isothiocyanate(smi)
        print(f"{name}: {result}, Reason: {reason}")