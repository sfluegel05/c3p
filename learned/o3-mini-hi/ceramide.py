"""
Classifies: CHEBI:17761 ceramide
"""
#!/usr/bin/env python3
"""
Classifies: Ceramide (N-acyl-sphingoid bases)

Ceramides are sphingoid base derivatives that contain an amide‐linked fatty acid.
The fatty acid is typically saturated or monounsaturated with chain lengths between 14 and 26 carbons.
Many sphingoid bases display a hydroxyl group (often on carbon 2, i.e. adjacent to the amide nitrogen).

This implementation:
  • Parses the molecule and adds explicit hydrogens (to improve -OH detection)
  • Looks for amide groups using a SMARTS pattern.
  • For each amide, from the carbonyl (C=O) side we “crawl” a contiguous, acyclic carbon chain
    requiring a chain length between 14 and 26 carbons.
  • From the nitrogen side, we check if one of the neighboring carbons (other than the carbonyl carbon)
    carries an -OH group.
    
Note: Due to the SMARTS "C(=O)N", the match will include three atoms (carbon, oxygen, nitrogen).
      We address that by extracting the first (carbonyl) and last (nitrogen) indices.
"""

from rdkit import Chem

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.

    It performs the following checks:
      1. The molecule must contain at least one amide group (C(=O)N).
      2. On the carbonyl side (fatty acyl side) of the amide, the chain must be a contiguous,
         acyclic chain (only following carbon–carbon bonds) with 14–26 carbon atoms.
      3. On the nitrogen side (sphingoid base), at least one carbon attached to the amide N must
         have a directly attached -OH group.
         
    Args:
       smiles (str): SMILES string for the molecule.

    Returns:
       bool: True if the molecule is classified as a ceramide, False otherwise.
       str: Explanation for the classification.
    """
    # Parse the SMILES string and add explicit hydrogens (better for detecting -OH groups)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # --- 1. Look for an amide group ---
    # The SMARTS "C(=O)N" matches three atoms: carbon, oxygen, nitrogen.
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide group (C(=O)N) found"
    
    # Helper: recursively determine the length of a contiguous, acyclic carbon chain.
    def linear_chain_length(atom_idx, coming_from, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Stop if the current atom is part of any ring.
        if atom.IsInRing():
            return 0
        max_length = 1
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == coming_from:
                continue
            # Only follow carbons that are not in a ring and haven’t been visited
            if nb.GetAtomicNum() == 6 and nb_idx not in visited:
                new_visited = visited.copy()
                new_visited.add(nb_idx)
                branch = 1 + linear_chain_length(nb_idx, atom_idx, new_visited)
                if branch > max_length:
                    max_length = branch
        return max_length

    # Helper: check for the sphingoid feature.
    # Look at neighbors of the amide nitrogen (except the carbonyl carbon) to see if any carbon has an -OH.
    def has_sphingoid_feature(n_atom, carbonyl_idx):
        for nb in n_atom.GetNeighbors():
            if nb.GetIdx() == carbonyl_idx:
                continue  # Skip the carbonyl carbon
            if nb.GetAtomicNum() == 6:  # Candidate carbon from the sphingoid base core
                # Check if any neighbor of this carbon is an oxygen with at least one hydrogen attached.
                for sub_nb in nb.GetNeighbors():
                    if sub_nb.GetAtomicNum() == 8:
                        bond = mol.GetBondBetweenAtoms(nb.GetIdx(), sub_nb.GetIdx())
                        if bond is None or bond.GetBondTypeAsDouble() != 1:
                            continue
                        # If the oxygen has at least one explicit hydrogen, count it as an -OH.
                        if sub_nb.GetTotalNumHs() > 0:
                            return True
        return False

    fatty_acyl_found = False
    sphingo_found = False
    reasons = []
    
    # Process each amide match.
    for match in amide_matches:
        # Because the SMARTS "C(=O)N" returns a tuple of three atoms,
        # we extract the carbonyl (first atom) and the nitrogen (last atom).
        carbonyl_idx = match[0]
        n_idx = match[-1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # --- 2. Identify the fatty acyl chain from the carbonyl side ---
        candidate_chain_idxs = []
        for nb in carbonyl_atom.GetNeighbors():
            if nb.GetIdx() == n_idx:
                continue  # Skip the amide nitrogen
            if nb.GetAtomicNum() != 6:
                continue  # Skip non-carbon atoms
            candidate_chain_idxs.append(nb.GetIdx())
        
        valid_chain_length = None
        for cand_idx in candidate_chain_idxs:
            cand_atom = mol.GetAtomWithIdx(cand_idx)
            if cand_atom.IsInRing():
                continue
            length = linear_chain_length(cand_idx, carbonyl_idx, {cand_idx})
            if 14 <= length <= 26:
                valid_chain_length = length
                fatty_acyl_found = True
                break
        if not fatty_acyl_found:
            continue  # Try next amide if no valid fatty acyl chain is found here
        
        # --- 3. Check for sphingoid side feature: an adjacent carbon with an -OH ---
        if has_sphingoid_feature(n_atom, carbonyl_idx):
            sphingo_found = True
            return True, ("Molecule is a ceramide: found an amide group with a fatty acyl chain of ~{} carbons "
                          "and a sphingoid base feature (N-linked carbon with an -OH).".format(valid_chain_length))
        else:
            reasons.append("Found an amide linkage with a valid fatty acyl chain (~{} carbons), "
                           "but the amide nitrogen is not attached to a carbon bearing an -OH.".format(valid_chain_length))
    
    if not fatty_acyl_found:
        return False, "No fatty acyl chain with 14–26 contiguous acyclic carbons (attached to an amide carbonyl) found"
    if not sphingo_found:
        return False, "Fatty acyl chain detected, but no suitable sphingoid base feature (N-adjacent carbon with -OH) was found"
    
    # Fallback (should not be reached)
    return False, "Structure does not meet ceramide criteria"

# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"  # e.g. N-eicosanoylsphinganine
# result, reason = is_ceramide(test_smiles)
# print(result, reason)