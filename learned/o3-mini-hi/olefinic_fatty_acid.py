"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid (free acid or as an acyl chain) that contains at least one C=C double bond.
A fatty acid is defined here as a linear (nearly unbranched) chain of at least 6 contiguous carbons attached
to a free acid group or as an acyl chain (e.g. esterified) that qualifies by having at least one carbon–carbon double bond.
We improve upon the previous DFS by enforcing linearity (no branching of the carbon chain).
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def dfs_linear_chain(mol, current_idx, previous_idx, current_length, found_double, min_length, visited):
    """
    Recursively traverse a strictly linear chain of sp3/sp2 carbons (ignoring ring atoms).
    At each step, if more than one candidate carbon to continue the chain is found, the chain is considered branched.
    A bond is checked for a double bond (C=C).
    Args:
        mol: RDKit molecule.
        current_idx: current atom idx (carbon).
        previous_idx: the atom idx from which we came.
        current_length: number of carbons in the chain so far.
        found_double: boolean flag indicating if a C=C has been encountered in the chain.
        min_length: minimum number of contiguous carbons required.
        visited: set of atom indices already visited.
    Returns:
        True if a linear chain of at least min_length (including current atom) is found and at least one double bond is encountered.
    """
    # If chain length (so far) is sufficient and we have seen a C=C, we accept.
    if current_length >= min_length and found_double:
        return True
    current_atom = mol.GetAtomWithIdx(current_idx)
    # Gather candidate neighbor carbons (ignoring the atom we came from and those in rings)
    candidates = []
    for nbr in current_atom.GetNeighbors():
        if nbr.GetAtomicNum() != 6:
            continue
        if nbr.GetIdx() == previous_idx:
            continue
        # Ignore if the neighbor is in a ring
        if nbr.IsInRing():
            continue
        candidates.append(nbr)
    # Enforce linearity: if more than one candidate to continue, then the chain is branched.
    if len(candidates) > 1:
        return False
    if len(candidates) == 0:
        # End of linear chain reached; check if the conditions are met.
        return False
    # Exactly one candidate: follow that bond.
    next_atom = candidates[0]
    # Get the bond between current and candidate to update found_double flag.
    bond = mol.GetBondBetweenAtoms(current_idx, next_atom.GetIdx())
    # Update flag if bond is a double bond.
    new_found_double = found_double or (bond and bond.GetBondType() == rdchem.BondType.DOUBLE)
    # Avoid cycles.
    if next_atom.GetIdx() in visited:
        return False
    visited.add(next_atom.GetIdx())
    result = dfs_linear_chain(mol, next_atom.GetIdx(), current_idx, current_length + 1, new_found_double, min_length, visited)
    visited.remove(next_atom.GetIdx())
    return result

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule qualifies as an olefinic fatty acid.
    The rules:
      1. Look for a free acid substructure (C(=O)[O;H1,O-]) or an acyl ester motif (C(=O)O[C]).
      2. For a free acid, the candidate “handle” is the alpha carbon neighbor of the carbonyl carbon.
         For an acyl ester, the handle is the carbon attached to the ester oxygen (in pattern C(=O)O[C]).
      3. From the handle, perform a DFS that enforces a linear (unbranched) chain of at least 6 contiguous carbons.
         During traversal, at least one C=C double bond must be encountered.
    Args:
         smiles: SMILES string of the molecule.
    Returns:
         (bool, str): True and a reason if classification is positive;
                      False and a reason message if not.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    MIN_CHAIN_LENGTH = 6  # at least 6 contiguous carbon atoms

    # Global quick check: molecule must include at least one carbon–carbon double bond.
    dbl_bond_smarts = Chem.MolFromSmarts("[#6]=[#6]")
    if not mol.HasSubstructMatch(dbl_bond_smarts):
        return False, "No carbon–carbon double bond (C=C) found in molecule"
    
    reasons = []
    
    # First: search for a free acid motif.
    # Free acid defined (here) as C(=O)[O;H1,O-].
    free_acid_smarts = "C(=O)[O;H1,O-]"
    free_acid_pattern = Chem.MolFromSmarts(free_acid_smarts)
    free_matches = mol.GetSubstructMatches(free_acid_pattern)
    if free_matches:
        for match in free_matches:
            # In the free acid SMARTS, match[0] is the carbonyl carbon.
            acid_carbon = mol.GetAtomWithIdx(match[0])
            # Look for alpha carbon neighbors (must be carbon and not in a ring).
            alpha_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6 and not nbr.IsInRing()]
            if not alpha_neighbors:
                reasons.append("Free acid group with no carbon neighbor found")
                continue
            for alpha in alpha_neighbors:
                # Initialize visited with acid carbon and the alpha candidate.
                visited = {acid_carbon.GetIdx(), alpha.GetIdx()}
                # For free acid the connecting bond is normally single; start DFS with current_length=1.
                if dfs_linear_chain(mol, alpha.GetIdx(), acid_carbon.GetIdx(), 1, False, MIN_CHAIN_LENGTH, visited):
                    return True, "Contains a fatty acyl chain (via free acid) with sufficient linear length and a C=C double bond."
            reasons.append("None of the free acid chains qualified (either too short, branched, or lacking a C=C double bond)")
    else:
        reasons.append("No free acid substructure (C(=O)[O;H1,O-]) found")
    
    # Second: search for acyl ester motifs.
    # Acyl ester motif: C(=O)O[C] where handle is the carbon attached to O.
    acyl_ester_smarts = "C(=O)O[C]"
    acyl_ester_pattern = Chem.MolFromSmarts(acyl_ester_smarts)
    ester_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    if ester_matches:
        for match in ester_matches:
            # For the ester motif, match[0] is carbonyl carbon, match[1] is oxygen, match[2] is the handle carbon.
            if len(match) < 3:
                continue
            handle_atom = mol.GetAtomWithIdx(match[2])
            # Start DFS from handle_atom.
            visited = {match[1], handle_atom.GetIdx()}
            if dfs_linear_chain(mol, handle_atom.GetIdx(), match[1], 1, False, MIN_CHAIN_LENGTH, visited):
                return True, "Contains a fatty acyl chain (via acyl ester) with sufficient linear length and a C=C double bond."
        reasons.append("Acyl ester chains did not qualify (either too short, branched, or lacking a C=C double bond)")
    else:
        reasons.append("No acyl ester substructure (C(=O)O[C]) found")
    
    # If no candidate chain qualifies, compile and return reasons.
    return False, " ; ".join(reasons)

# Example usage:
if __name__ == "__main__":
    test_smiles = [
        # True positives:
        "O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O",  # Avenoleic acid
        "O=C(CCC/C=C\\C/C=C\\CC(/C=C/C(C(CCCCC)O)O)O)O",  # 11,14,15-trihydroxy-(5Z,8Z,12E)-icosatrienoic acid
        "CC\\C=C/C\\C=C/CCC\\C=C\\C=C\\C=C/CCCC(O)=O",  # (5Z,7E,9E,14Z,17Z)-icosapentaenoic acid
        "[H]C(\\C=C/CCCCCCCC(O)=O)=C1O[C@H]1C\\C=C/CC",  # (9Z,13S,15Z)-12,13-epoxyoctadeca-9,11,15-trienoic acid
        "OC(CCCCC)\\C=C\\C=C\\CCCCCCCC(O)=O",  # alpha-Artemisic acid
        "CCCCCCCCC\\C=C/CCCCCCCC(O)=O",  # gadoleic acid (shortened simulated example)
        # False positives (should be classified as False):
        "O(C(C[N+](C)(C)C)CC([O-])=O)C(=O)CC/C=C/C/C=C/C/C=C/C\\C=C\\CC(O)/C=C/C=C/CC",  # Complex acyl chain in a non-fatty acid moiety
        # False negatives (should be classified as True but were missed previously – if the chain is isolated and meets the criteria):
        "OC[C@H](COP(OCC[N+](C)(C)C)(=O)[O-])OC(=O)CCCCCCC/C=C\\CCCCCC",  # 2-[(9Z)-hexadecenoyl]-sn-glycero-3-phosphocholine
    ]
    
    for sm in test_smiles:
        result, reason = is_olefinic_fatty_acid(sm)
        print("SMILES:", sm)
        print("Result:", result)
        print("Reason:", reason)
        print("-" * 60)