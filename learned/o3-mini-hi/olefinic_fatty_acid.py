"""
Classifies: CHEBI:53339 olefinic fatty acid
"""
"""
Classifies: Olefinic Fatty Acid
Definition: Any fatty acid (free or as an acyl chain) that 
contains at least one C=C double bond.
A fatty acid is defined as a linear, unbranched (or nearly unbranched) chain 
of at least 6 contiguous carbon atoms that has at least one carbon–carbon double bond.
This implementation uses a DFS of the candidate chain (ignoring atoms in rings)
to find any simple (non‐cyclic) path that is long enough and contains a double bond.
"""
from rdkit import Chem
from rdkit.Chem import rdchem

# Depth-first search to try to find a linear carbon chain (path) of minimum length with a C=C.
def dfs_chain(mol, current_atom, prev_idx, current_length, has_double, min_length, visited):
    # If we already have a long enough chain with a double bond, we can return True.
    if current_length >= min_length and has_double:
        return True
    # For each neighbor that is carbon, not coming back, and not in a ring:
    for nbr in current_atom.GetNeighbors():
        # We want only carbon atoms and avoid revisiting the previous atom (and already visited atoms).
        if nbr.GetAtomicNum() != 6 or nbr.GetIdx() == prev_idx or nbr.GetIdx() in visited:
            continue
        # Also ignore atoms that are part of a ring
        if nbr.IsInRing():
            continue
        bond = mol.GetBondBetweenAtoms(current_atom.GetIdx(), nbr.GetIdx())
        # If the bond is a double bond, mark it.
        new_double = has_double or (bond and bond.GetBondType() == rdchem.BondType.DOUBLE)
        visited.add(nbr.GetIdx())
        if dfs_chain(mol, nbr, current_atom.GetIdx(), current_length + 1, new_double, min_length, visited):
            return True
        visited.remove(nbr.GetIdx())
    return False

def is_olefinic_fatty_acid(smiles: str):
    """
    Determines if a molecule is an olefinic fatty acid.

    The rules:
      1. Identify a free acid or acyl ester substructure.
      2. For each, the carbon atom attached to the acid/ester group is considered the handle.
      3. From the handle, search (via DFS) among non-ring, carbon–carbon bonds for a simple path
         of at least 6 contiguous carbons (including the handle) that also contains at least one double bond.
      4. If such a path exists, then the molecule qualifies.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as an olefinic fatty acid, False otherwise.
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    MIN_CHAIN_LENGTH = 6

    # Quick global check: must be at least one C=C in the whole molecule.
    dbl_bond_global = Chem.MolFromSmarts("[#6]=[#6]")
    if not mol.HasSubstructMatch(dbl_bond_global):
        return False, "No carbon–carbon double bond (C=C) found in molecule"
    
    reasons = []  # To collect error messages along each attempted route

    # First, try free acid motif.
    # Free acid defined as C(=O)[O;H1,O-]
    free_acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    free_matches = mol.GetSubstructMatches(free_acid_pattern)
    if free_matches:
        for match in free_matches:
            # In the free acid SMARTS, match[0] is the carbonyl carbon.
            acid_c = mol.GetAtomWithIdx(match[0])
            # Find carbon neighbors of the acid carbon (the "alpha" carbon(s)).
            alpha_neighbors = [nbr for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
            if not alpha_neighbors:
                reasons.append("Free acid group with no carbon neighbor found")
                continue
            # For each candidate alpha carbon, search for a fatty acyl chain.
            for alpha in alpha_neighbors:
                # Do not require the alpha to be strictly terminal.
                visited = set([acid_c.GetIdx(), alpha.GetIdx()])
                if dfs_chain(mol, alpha, acid_c.GetIdx(), 1, False, MIN_CHAIN_LENGTH, visited):
                    return True, "Contains a fatty acyl chain (via free acid) with sufficient length and a C=C double bond."
            reasons.append("None of the free acid chains qualified (either too short or lacking a C=C double bond)")
    
    # Next, try acyl ester motif.
    # Acyl ester motif: C(=O)O[C] . Here the handle is the carbon attached to the ester oxygen.
    acyl_ester_pattern = Chem.MolFromSmarts("C(=O)O[C]")
    ester_matches = mol.GetSubstructMatches(acyl_ester_pattern)
    if ester_matches:
        for match in ester_matches:
            # match[0]: carbonyl carbon; match[1]: oxygen; match[2]: handle carbon.
            if len(match) < 3:
                continue
            handle_atom = mol.GetAtomWithIdx(match[2])
            # In phospholipids the handle might be slightly branched.
            visited = set([match[1], handle_atom.GetIdx()])
            if dfs_chain(mol, handle_atom, match[1], 1, False, MIN_CHAIN_LENGTH, visited):
                return True, "Contains a fatty acyl chain (via acyl ester) with sufficient length and a C=C double bond."
            reasons.append("Acyl ester chain starting at handle did not qualify (either too short or lacking C=C)")
    
    # If nothing qualifies, compile the reasons.
    if reasons:
        return False, " ; ".join(reasons)
    return False, "No free acid or acyl ester substructure with a qualifying fatty acyl chain found."

# Example usage: (for testing; remove or comment out during integration)
if __name__ == "__main__":
    test_smiles_list = [
        # True positives
        "O[C@H](CCC)C/C=C\\C/C=C\\CCCCCCCC(O)=O",  # Avenoleic acid
        "O=C(CCC/C=C\\C/C=C\\CC(/C=C/C(C(CCCCC)O)O)O)O",  # 11,14,15-trihydroxy-(5Z,8Z,12E)-icosatrienoic acid
        "CC\\C=C/C\\C=C/CCC\\C=C\\C=C\\C=C/CCCC(O)=O",  # (5Z,7E,9E,14Z,17Z)-icosapentaenoic acid
        "[H]C(\\C=C/CCCCCCCC(O)=O)=C1O[C@H]1C\\C=C/CC",  # (9Z,13S,15Z)-12,13-epoxyoctadeca-9,11,15-trienoic acid
        "OC(CCCCC)\\C=C\\C=C\\CCCCCCCC(O)=O",  # alpha-Artemisic acid
        # False positives (should be False)
        "O=C(O)/C=C/C=C/[C@@H]([C@@H]1O[C@]2(O[C@@H]([C@@H](CC)[C@@H](C2)O)C)[C@@H](C)[C@@H]([C@@H]1C)O)C",  # Pteridic acid F
        # False negatives (should be True, but were missed previously)
        "OC[C@H](COP(OCC[N+](C)(C)C)(=O)[O-])OC(=O)CCCCCCC/C=C\\CCCCCC",  # 2-[(9Z)-hexadecenoyl]-sn-glycero-3-phosphocholine
    ]
    
    for s in test_smiles_list:
        res, msg = is_olefinic_fatty_acid(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n{'-'*60}")