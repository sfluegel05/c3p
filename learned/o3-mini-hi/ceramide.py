"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: Ceramide (N-acyl-sphingoid bases)

Ceramides are sphingoid base derivatives that contain an amide‐linked fatty acid.
The fatty acid is typically saturated or monounsaturated with chain lengths between 14 and 26 carbons.
Many sphingoid bases display a hydroxyl group (often on carbon 2, i.e. adjacent to the amide nitrogen).

This implementation:
  • Parses the molecule (with added explicit hydrogens for -OH detection)
  • Looks for amide groups using a SMARTS pattern.
  • For each amide, from the carbonyl (C=O) side we try to “crawl” a contiguous, acyclic carbon chain;
    we require a chain length between 14 and 26 carbons.
  • From the nitrogen side we check for a sphingoid base feature by inspecting the carbon(s)
    attached to the amide nitrogen (except the carbonyl) to see if one has a direct –OH substituent.
  
Note: This heuristic may mis‐classify unusual structures.
"""

from rdkit import Chem

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide, based on its SMILES string.

    It performs the following checks:
      1. The molecule must contain at least one amide group (C(=O)N).
      2. On the carbonyl side (fatty acyl side) of the amide, there should be a contiguous,
         acyclic carbon chain of 14–26 atoms. (Only following carbon–carbon bonds that are not in rings.)
      3. On the nitrogen side (sphingoid base), the amide N should be attached to at least one carbon
         that carries an -OH group (as a proxy for the common hydroxyl on carbon 2).
         
    Args:
       smiles (str): SMILES string for the molecule

    Returns:
       bool: True if the molecule is classified as a ceramide, False otherwise.
       str: Explanation for the classification.
    """
    # Parse SMILES and add explicit hydrogens (improves -OH detection)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # --- 1. Look for an amide group ---
    # This SMARTS will match a carbonyl carbon directly joined to a nitrogen.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide group (C(=O)N) found"
    
    # Helper: recursively determine the length of a contiguous, acyclic carbon chain.
    def linear_chain_length(atom_idx, coming_from, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        # do not travel further if the atom is in a ring
        if atom.IsInRing():
            return 0
        max_length = 1
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == coming_from:
                continue
            # only follow carbons that are not in a ring and haven’t been visited
            if nb.GetAtomicNum() == 6 and nb_idx not in visited:
                new_visited = visited.copy()
                new_visited.add(nb_idx)
                branch = 1 + linear_chain_length(nb_idx, atom_idx, new_visited)
                if branch > max_length:
                    max_length = branch
        return max_length

    # Helper: check for the sphingoid base feature.
    # Here we examine neighbors of the amide N (except the carbonyl atom) to see if any neighboring carbon
    # carries at least one –OH group attached.
    def has_sphingoid_feature(n_atom, carbonyl_idx):
        for nb in n_atom.GetNeighbors():
            if nb.GetIdx() == carbonyl_idx:
                continue  # skip the carbonyl carbon
            if nb.GetAtomicNum() == 6:  # candidate carbon from the sphingoid base core
                # check its neighbors for an -OH (oxygen with a single bond and at least one hydrogen)
                for sub_nb in nb.GetNeighbors():
                    if sub_nb.GetAtomicNum() == 8:
                        # Check bond type between nb and sub_nb (should be a single bond)
                        bond = mol.GetBondBetweenAtoms(nb.GetIdx(), sub_nb.GetIdx())
                        if bond is None or bond.GetBondTypeAsDouble() != 1:
                            continue
                        # If the oxygen has at least one explicit hydrogen, count it as -OH.
                        if sub_nb.GetTotalNumHs() > 0:
                            return True
        return False

    fatty_acyl_found = False
    sphingo_found = False
    reasons = []
    
    # Process each amide match: each match is a tuple of (carbonyl_idx, nitrogen_idx)
    for match in amide_matches:
        carbonyl_idx, n_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # --- 2. Identify the fatty acyl chain from the carbonyl side. ---
        # Look at all neighbors of the carbonyl except the amide nitrogen.
        candidate_chain_idxs = []
        for nb in carbonyl_atom.GetNeighbors():
            if nb.GetIdx() == n_idx:
                continue
            # Skip oxygen atoms (e.g., the carbonyl O or other oxygens)
            if nb.GetAtomicNum() != 6:
                continue
            candidate_chain_idxs.append(nb.GetIdx())
        
        valid_chain_length = None
        for cand_idx in candidate_chain_idxs:
            cand_atom = mol.GetAtomWithIdx(cand_idx)
            # Do not follow if candidate carbon is in a ring
            if cand_atom.IsInRing():
                continue
            length = linear_chain_length(cand_idx, carbonyl_idx, {cand_idx})
            if 14 <= length <= 26:
                valid_chain_length = length
                fatty_acyl_found = True
                break
        if not fatty_acyl_found:
            continue  # try the next amide if no valid fatty acyl chain found for this match
        
        # --- 3. Check for sphingoid feature on the nitrogen side. ---
        if has_sphingoid_feature(n_atom, carbonyl_idx):
            sphingo_found = True
            return True, ("Molecule is a ceramide: found an amide group with an acyl chain of ~{} carbons "
                          "and a sphingoid base feature (N-linked carbon with an -OH).".format(valid_chain_length))
        else:
            reasons.append("Found an amide linkage with a valid fatty acyl chain (~{} carbons), "
                           "but the amide nitrogen is not attached to a carbon bearing -OH.".format(valid_chain_length))
    
    # Give a reason based on what was (or was not) found.
    if not fatty_acyl_found:
        return False, "No fatty acyl chain with 14–26 contiguous acyclic carbons (attached to an amide carbonyl) found"
    if not sphingo_found:
        return False, "Fatty acyl chain detected, but no suitable sphingoid base feature (N-adjacent carbon with -OH) was found"
    
    # Fallback reason (should not reach here)
    return False, "Structure does not meet ceramide criteria"

# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"  # e.g. N-eicosanoylsphinganine
# result, reason = is_ceramide(test_smiles)
# print(result, reason)