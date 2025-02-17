"""
Classifies: CHEBI:17761 ceramide
"""
#!/usr/bin/env python3
"""
Classifies: Ceramide (N-acyl-sphingoid bases)

Ceramides are sphingoid base derivatives characterized by an amide‐linked fatty acid.
The fatty acid is typically saturated or monounsaturated and has a chain length of 14–26 carbons.
A hydroxyl group on the sphingoid base (often near the amide nitrogen, e.g. on C2) is common.
This implementation:
  • Parses the molecule and adds explicit hydrogens (to improve –OH detection)
  • Searches for amide groups (using SMARTS “C(=O)N”)
  • For each amide, examines the carbonyl side for a contiguous, acyclic carbon chain of 14–26 atoms.
  • From the amide nitrogen side, performs a breadth-first search (up to 3 bonds) for an oxygen
    atom with an –OH moiety.
If both features are found on any of the amide linkages, the molecule is classified as a ceramide.
"""

from rdkit import Chem

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide based on its SMILES string.
    
    It performs the following checks for each amide linkage (C(=O)N):
      1. On the carbonyl side (fatty acyl side), there must be a contiguous, acyclic carbon chain
         of 14–26 carbons.
      2. On the nitrogen side (sphingoid base), within up to 3 bonds (excluding the carbonyl),
         an oxygen atom with at least one hydrogen (likely an –OH group) must be found.
    
    Args:
       smiles (str): SMILES string representing the molecule.
    
    Returns:
       bool: True if classified as a ceramide, False otherwise.
       str: Explanation for the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens to aid in detecting –OH groups
    mol = Chem.AddHs(mol)
    
    # 1. Look for an amide linkage (SMARTS: C(=O)N)
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide group (C(=O)N) found"
    
    # Helper: recursively determine the maximum length of a contiguous, acyclic carbon chain
    def linear_chain_length(atom_idx, coming_from, visited):
        atom = mol.GetAtomWithIdx(atom_idx)
        # Do not follow through rings
        if atom.IsInRing():
            return 0
        max_length = 1
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == coming_from:
                continue
            if nb.GetAtomicNum() == 6 and nb_idx not in visited:
                new_visited = visited.copy()
                new_visited.add(nb_idx)
                branch_length = 1 + linear_chain_length(nb_idx, atom_idx, new_visited)
                if branch_length > max_length:
                    max_length = branch_length
        return max_length
    
    # Helper: from the amide nitrogen, perform a BFS (up to depth 3) excluding the carbonyl atom,
    # to search for an oxygen atom that appears to be in an –OH group (i.e. has at least one hydrogen).
    def has_sphingoid_feature(n_atom, carbonyl_idx):
        from collections import deque
        visited = set()
        # Queue contains tuples: (atom, depth)
        queue = deque()
        # Start with neighbors of the amide nitrogen (except the carbonyl)
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetIdx() == carbonyl_idx:
                continue
            queue.append((neighbor, 1))
            visited.add(neighbor.GetIdx())
        max_depth = 3
        while queue:
            current, depth = queue.popleft()
            # Check if current atom is oxygen with an –OH (i.e. at least one explicit hydrogen)
            if current.GetAtomicNum() == 8 and current.GetTotalNumHs() > 0:
                return True
            # Otherwise, if we have not reached max depth, add neighbors (except we do not go back to n_atom)
            if depth < max_depth:
                for nb in current.GetNeighbors():
                    nb_idx = nb.GetIdx()
                    if nb_idx in visited:
                        continue
                    # Avoid going back to the original amide nitrogen or to the carbonyl atom
                    if nb.GetIdx() == n_atom.GetIdx() or nb.GetIdx() == carbonyl_idx:
                        continue
                    visited.add(nb_idx)
                    queue.append((nb, depth + 1))
        return False

    # Process each amide match
    for match in amide_matches:
        # match is a tuple of atom indices corresponding to (carbonyl carbon, oxygen, amide N)
        # We only need the carbonyl carbon (first atom) and the amide nitrogen (last atom)
        carbonyl_idx = match[0]
        n_idx = match[-1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # --- 2. Check the fatty acyl chain: examine neighbors of the carbonyl (excluding the amide N)
        fatty_acyl_found = False
        valid_chain_length = None
        for nb in carbonyl_atom.GetNeighbors():
            if nb.GetIdx() == n_idx:  # skip the amide nitrogen
                continue
            if nb.GetAtomicNum() != 6:  # must be carbon
                continue
            # Only consider chains that are not in rings.
            if nb.IsInRing():
                continue
            # Use our DFS helper to count the number of contiguous acyclic carbons:
            chain_len = linear_chain_length(nb.GetIdx(), carbonyl_idx, {nb.GetIdx()})
            if 14 <= chain_len <= 26:
                fatty_acyl_found = True
                valid_chain_length = chain_len
                break
        
        if not fatty_acyl_found:
            # Try next amide match if available
            continue
        
        # --- 3. Check the sphingoid base side: from the nitrogen, search for an –OH within up to 3 bonds.
        if has_sphingoid_feature(n_atom, carbonyl_idx):
            message = ("Molecule is classified as a ceramide: Amide linkage found with a fatty acyl chain "
                       "of ~{} contiguous carbons and a nearby sphingoid base feature (–OH within three bonds "
                       "of the amide nitrogen).".format(valid_chain_length))
            return True, message
        # If an amide match had a valid fatty acyl chain but the sphingoid base feature is missing,
        # continue to see if another amide group meets both criteria.
    
    # If no amide linkage satisfies both criteria, return an explanation.
    # Determine which feature was never found.
    if not any(linear_chain_length(nb.GetIdx(), match[0], {nb.GetIdx()}) >= 14 and 
               linear_chain_length(nb.GetIdx(), match[0], {nb.GetIdx()}) <= 26 
               for match in amide_matches 
               for nb in mol.GetAtomWithIdx(match[0]).GetNeighbors() if nb.GetIdx() != match[-1] and nb.GetAtomicNum() == 6):
        return False, "No fatty acyl chain with 14–26 contiguous acyclic carbons (attached to an amide carbonyl) found"
    else:
        return False, "Fatty acyl chain found, but no suitable sphingoid base feature (–OH near the amide N) was detected"

# Example usage:
# Uncomment the following lines to test with one of the ceramide SMILES examples:
# test_smiles = "CCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"  # N-eicosanoylsphinganine
# result, reason = is_ceramide(test_smiles)
# print(result, reason)