"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate (R-N=C=S, an organosulfur compound)

An isothiocyanate is defined as an organosulfur compound with the general formula R-N=C=S.
This improved function first looks for the characteristic SMARTS pattern (in both directions)
and then enforces stricter connectivity rules by checking:
  - Both bonds (N–C and C–S) are double bonds.
  - The central C atom has exactly 2 heavy neighbors (only N and S).
  - The terminal S atom has exactly 1 heavy neighbor (the central C).
  - The N atom has exactly 2 heavy neighbors (the central C and one substituent R).
  - Additionally, if the substituent R (bonded to N) is aromatic and directly attached to a carbonyl group,
    the candidate is rejected.
If a candidate match passes these tests, it is assumed to be a terminal isothiocyanate group.
Otherwise the function returns False.
"""

from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule belongs to the isothiocyanate class based on its SMILES string.
    An isothiocyanate is defined as an organosulfur compound with the general formula R-N=C=S.
    This function looks for the SMARTS pattern "N=C=S" (or its reverse "S=C=N") and then enforces
    stricter connectivity checks so that:
      - the central C is only bonded to N and S,
      - the terminal S is bonded solely to that C, and
      - the N is bonded to exactly one substituent group (R) besides the C.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as an isothiocyanate, False otherwise.
        str: Reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define two SMARTS patterns for the isothiocyanate group:
    # pattern1: in the order N=C=S
    pattern1 = Chem.MolFromSmarts("N=C=S")
    # pattern2: in reverse order S=C=N (will be re-ordered later)
    pattern2 = Chem.MolFromSmarts("S=C=N")
    
    candidate_matches = []
    if mol.HasSubstructMatch(pattern1):
        candidate_matches.extend(mol.GetSubstructMatches(pattern1))
    if mol.HasSubstructMatch(pattern2):
        candidate_matches.extend(mol.GetSubstructMatches(pattern2))
    
    for match in candidate_matches:
        # Depending on which pattern was detected the atom order may differ.
        # Default assume match from pattern1 gives (N, C, S)
        atom0 = mol.GetAtomWithIdx(match[0])
        atom1 = mol.GetAtomWithIdx(match[1])
        atom2 = mol.GetAtomWithIdx(match[2])
        
        n_idx, c_idx, s_idx = match[0], match[1], match[2]
        # If first atom is S (atomic number 16), then assume the match is from the reverse pattern S=C=N.
        if atom0.GetAtomicNum() == 16:
            # Reorder the match to (N, C, S)
            n_idx, c_idx, s_idx = match[2], match[1], match[0]
            
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        s_atom = mol.GetAtomWithIdx(s_idx)
        
        # Check that the bonds between N-C and C-S exist and are double bonds.
        bond_nc = mol.GetBondBetweenAtoms(n_idx, c_idx)
        bond_cs = mol.GetBondBetweenAtoms(c_idx, s_idx)
        if bond_nc is None or bond_cs is None:
            continue
        if bond_nc.GetBondType() != Chem.rdchem.BondType.DOUBLE or bond_cs.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            continue
        
        # Enforce terminal connectivity using heavy-atom counts.
        # For the central C: only heavy neighbors should be n_atom and s_atom.
        heavy_neighbors_c = [nbr.GetIdx() for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors_c) != 2:
            continue
        
        # The terminal S should have only the central C as a heavy neighbor.
        heavy_neighbors_s = [nbr.GetIdx() for nbr in s_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors_s) != 1:
            continue
        
        # The N atom should have exactly two heavy neighbors: one is c_atom and one is the R substituent.
        heavy_neighbors_n = [nbr.GetIdx() for nbr in n_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors_n) != 2:
            continue
        
        # (Optional additional check)
        # Let r_idx be the neighbor of N that is not the central carbon.
        r_idx = None
        for idx in heavy_neighbors_n:
            if idx != c_idx:
                r_idx = idx
                break
        if r_idx is None:
            continue
        r_atom = mol.GetAtomWithIdx(r_idx)
        
        # If the R substituent is aromatic, we check if it is directly attached to a carbonyl,
        # which might indicate that the isothiocyanate group is part of a larger conjugated system.
        if r_atom.GetIsAromatic():
            attached_to_carbonyl = False
            for nbr in r_atom.GetNeighbors():
                # Skip if neighbor is N (the one already used)
                if nbr.GetIdx() == n_idx:
                    continue
                # Look for a carbon (atomic number 6) that has a double bond to oxygen.
                if nbr.GetAtomicNum() == 6:
                    for bond in nbr.GetBonds():
                        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                            other = bond.GetOtherAtom(nbr)
                            if other.GetAtomicNum() == 8:
                                attached_to_carbonyl = True
                                break
                if attached_to_carbonyl:
                    break
            if attached_to_carbonyl:
                continue  # Reject candidate
        
        # Also enforce that the degrees (number of bonds) of the three atoms are exactly as expected.
        # For an ideal terminal isothiocyanate we expect:
        #   N: degree 2 (bond to R and C)
        #   C: degree 2 (bond to N and S)
        #   S: degree 1 (bond to C)
        if n_atom.GetDegree() != 2 or c_atom.GetDegree() != 2 or s_atom.GetDegree() != 1:
            continue
        
        # If all tests are passed, then we have found a terminal isothiocyanate group.
        return True, "Contains terminal isothiocyanate group (R-N=C=S)"
    
    return False, "Does not contain a suitable terminal isothiocyanate group"


# Example usage and quick tests
if __name__ == "__main__":
    # List a few example SMILES strings (true and false candidates)
    test_smiles = {
        "sulforaphane": "N(=C=S)CCCCS(=O)C",
        "phenyl isothiocyanate": "S=C=Nc1ccccc1",
        "Camelinin (false positive candidate)": "S(=O)(CCCCCCCCCCN=C=S)C",
        "4-Isothiocyanatobenzoate (false positive candidate)": "[O-]C(=O)c1ccc(cc1)N=C=S"
    }
    
    for name, smi in test_smiles.items():
        result, reason = is_isothiocyanate(smi)
        print(f"{name}: {result}, Reason: {reason}")