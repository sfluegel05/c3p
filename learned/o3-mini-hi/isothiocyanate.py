"""
Classifies: CHEBI:52221 isothiocyanate
"""
"""
Classifies: Isothiocyanate (R-N=C=S, an organosulfur compound)
An isothiocyanate is defined as an organosulfur compound with the general formula R-N=C=S.
This improved function not only searches for the characteristic SMARTS pattern(s)
but also checks the atomic connectivity, bond orders, and terminality of the group.
"""
from rdkit import Chem

def is_isothiocyanate(smiles: str):
    """
    Determines if a molecule belongs to the isothiocyanate class based on its SMILES string.
    An isothiocyanate is defined as an organosulfur compound with the general formula R-N=C=S.
    This function first searches for the SMARTS for "N=C=S" (and its reverse "S=C=N")
    and then inspects for a terminal (non-embedded) isothiocyanate group by checking:
      - The central carbon (C) has exactly two heavy neighbors (N and S only).
      - The sulfur (S) has only one heavy neighbor (the C).
      - Both the N-C and C-S bonds are double bonds.
    
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
    # pattern2: sometimes written in reverse as S=C=N
    pattern2 = Chem.MolFromSmarts("S=C=N")
    
    # Gather all candidate matches. Each match is a tuple of atom indices.
    candidate_matches = []
    if mol.HasSubstructMatch(pattern1):
        candidate_matches.extend(mol.GetSubstructMatches(pattern1))
    if mol.HasSubstructMatch(pattern2):
        candidate_matches.extend(mol.GetSubstructMatches(pattern2))
    
    # Process each candidate match.
    # We want to have the match in a consistent order as (n_idx, c_idx, s_idx)
    for match in candidate_matches:
        # Determine ordering based on which SMARTS matched.
        # For a match from pattern1 (N=C=S), we expect match = (n, c, s)
        # For pattern2 (S=C=N), the order is (s, c, n) so we reorder.
        atom0 = mol.GetAtomWithIdx(match[0])
        atom1 = mol.GetAtomWithIdx(match[1])
        atom2 = mol.GetAtomWithIdx(match[2])
        
        # By default assume ordering from pattern1:
        n_idx, c_idx, s_idx = match[0], match[1], match[2]
        # If the first atom is sulfur (atomic number 16) then we assume pattern2 order.
        if atom0.GetAtomicNum() == 16:
            # Then atom order in the match is (S, C, N). We reverse to (N, C, S)
            n_idx, c_idx, s_idx = match[2], match[1], match[0]
        
        n_atom = mol.GetAtomWithIdx(n_idx)
        c_atom = mol.GetAtomWithIdx(c_idx)
        s_atom = mol.GetAtomWithIdx(s_idx)
        
        # Get the bonds between N-C and C-S.
        bond_nc = mol.GetBondBetweenAtoms(n_idx, c_idx)
        bond_cs = mol.GetBondBetweenAtoms(c_idx, s_idx)
        if bond_nc is None or bond_cs is None:
            continue  # Not a proper match
        
        # Check that both bonds are double bonds
        if bond_nc.GetBondType() != Chem.rdchem.BondType.DOUBLE or bond_cs.GetBondType() != Chem.rdchem.BondType.DOUBLE:
            continue
        
        # Check that the isothiocyanate group is "terminal":
        # 1) The central carbon should have exactly 2 heavy-atom neighbors (only n and s).
        heavy_neighbors_c = [nbr.GetIdx() for nbr in c_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors_c) != 2:
            continue
        
        # 2) The sulfur atom should have only 1 heavy neighbor (the carbon).
        heavy_neighbors_s = [nbr.GetIdx() for nbr in s_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors_s) != 1:
            continue
        
        # 3) (Optional) Check that the nitrogen atom has at least one substituent (the R group)
        #    apart from the isothiocyanate carbon.
        heavy_neighbors_n = [nbr for nbr in n_atom.GetNeighbors() if nbr.GetIdx() != c_idx and nbr.GetAtomicNum() > 1]
        if len(heavy_neighbors_n) < 1:
            continue

        # This candidate passed our connectivity & bond order checks.
        return True, "Contains terminal isothiocyanate group (R-N=C=S)"
    
    return False, "Does not contain a suitable terminal isothiocyanate group"
    
# Example usage:
if __name__ == "__main__":
    # A few test SMILES for true positives and false positives
    test_smiles = {
        "sulforaphane": "N(=C=S)CCCCS(=O)C",
        "phenyl isothiocyanate": "S=C=Nc1ccccc1",
        "Camelinin (false positive candidate)": "S(=O)(CCCCCCCCCCN=C=S)C",
        "Hapalindole M (false positive candidate)": "S=C=N[C@H]1[C@@](C=C)(CC[C@H]2[C@@H]1C=3C4=C(C=CC=C4NC3)C2(C)C)C"
    }
    for name, smi in test_smiles.items():
        result, reason = is_isothiocyanate(smi)
        print(f"{name}: {result}, Reason: {reason}")