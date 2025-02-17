"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: Ceramide (N-acyl-sphingoid bases)

Ceramides are sphingoid base derivatives that contain an amide‐linked fatty acid.
The fatty acid is typically saturated or monounsaturated with chain lengths between 14 and 26 carbons.
Many sphingoid bases display a hydroxyl group (often on carbon 2, i.e. adjacent to the amide nitrogen).
This improved implementation:
  • Adds explicit hydrogens so that hydroxyl groups are detected reliably.
  • Checks for an amide group with a fatty acyl chain (14–26 contiguous non‐ring carbons)
    on the carbonyl (acyl) side.
  • Checks that the amide nitrogen (base side) is bonded to a carbon that directly bears an –OH group.
Note: This is heuristic and may mis‐classify unusual structures.
"""
from rdkit import Chem

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide, based on its SMILES string.

    It performs the following checks:
      1. The molecule must contain an amide group (C(=O)N).
      2. On the carbonyl side of an amide, there should be a contiguous, acyclic carbon chain of 14–26 atoms.
         (This chain represents the fatty acyl portion.)
      3. On the nitrogen side of that same amide, the nitrogen should be attached (directly) to a carbon that carries
         a hydroxyl group (-OH), which is a common feature of the sphingoid base.

    Args:
       smiles (str): SMILES string for the molecule

    Returns:
       bool: True if the molecule is classified as a ceramide, False otherwise.
       str: Explanation for the classification.
    """
    # Parse SMILES and add explicit hydrogens for reliable -OH detection.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # --- 1. Look for an amide group ---
    # Use a simple SMARTS for an amide. The match returns a tuple (carbonyl carbon, amide nitrogen)
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide group (C(=O)N) found"
    
    # Helper function: recursively determine the length of a contiguous, acyclic carbon chain.
    def linear_chain_length(mol, current_idx, coming_from_idx, visited):
        """
        Depth-first search (DFS) to return the maximum length of a contiguous, acyclic carbon chain.
        Only follows carbon atoms that are not part of a ring.
        """
        current_atom = mol.GetAtomWithIdx(current_idx)
        if current_atom.IsInRing():
            return 0
        max_length = 1
        for nb in current_atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == coming_from_idx:
                continue
            if nb.GetAtomicNum() != 6:
                continue
            if nb_idx in visited:
                continue
            new_visited = visited.copy()
            new_visited.add(nb_idx)
            branch_length = 1 + linear_chain_length(mol, nb_idx, current_idx, new_visited)
            if branch_length > max_length:
                max_length = branch_length
        return max_length

    # Helper function: check for a sphingoid base feature.
    def has_sphingoid_feature(n_atom, carbonyl_idx):
        """
        Checks that the amide nitrogen (n_atom) is attached to at least one carbon (other than the carbonyl carbon)
        which directly bears a hydroxyl (-OH) group.
        We assume that the molecule has explicit hydrogens so that an -OH will have at least one hydrogen.
        """
        for nb in n_atom.GetNeighbors():
            if nb.GetIdx() == carbonyl_idx:
                continue  # skip the carbonyl carbon
            if nb.GetAtomicNum() == 6:  # neighbor is carbon
                # Look among the neighbors of this carbon for an oxygen that is -OH.
                for sub_nb in nb.GetNeighbors():
                    # We want an oxygen with at least one hydrogen attached
                    if sub_nb.GetAtomicNum() == 8 and sub_nb.GetTotalNumHs() > 0:
                        return True
        return False

    fatty_acyl_found = False
    sphingo_found = False
    collected_reasons = []
    
    # --- 2 & 3. For each amide, check the acyl chain (from the carbonyl side) and then the sphingoid base feature.
    for match in amide_matches:
        carbonyl_idx = match[0]
        n_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # For the acyl (fatty acid) side, consider neighbors of the carbonyl carbon except the amide nitrogen
        candidate_chain_idxs = []
        for nb in carbonyl_atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx == n_idx:  # skip amide nitrogen
                continue
            # Skip typical carbonyl oxygen or any oxygen
            if nb.GetAtomicNum() == 8:
                continue
            if nb.GetAtomicNum() == 6:
                candidate_chain_idxs.append(nb_idx)
        
        # Check each candidate chain to see if its contiguous length is within 14 to 26 carbons.
        for cand in candidate_chain_idxs:
            if mol.GetAtomWithIdx(cand).IsInRing():
                continue
            length = linear_chain_length(mol, cand, carbonyl_idx, {cand})
            if 14 <= length <= 26:
                fatty_acyl_found = True
                break  # Found a valid fatty acyl chain for this amide link.
        if not fatty_acyl_found:
            continue  # Try the next amide match
        
        # Next, check if the amide nitrogen shows the sphingoid feature: i.e. attached carbon with a hydroxyl.
        if has_sphingoid_feature(n_atom, carbonyl_idx):
            sphingo_found = True
            return True, ("Molecule is a ceramide: it contains an amide-linked fatty acyl chain (~{} carbons) "
                          "and a sphingoid base feature (the amide N is attached to a carbon bearing an -OH)."
                          .format(length))
        else:
            collected_reasons.append("Amide group found with valid fatty acyl chain but no sphingoid base feature "
                                     "(N not connected to a carbon with -OH) for this linkage.")
    
    if not fatty_acyl_found:
        return False, "No fatty acyl chain with 14–26 contiguous carbons (attached to an amide carbonyl) found"
    if not sphingo_found:
        return False, "Fatty acyl chain detected, but no suitable sphingoid base feature (N-adjacent -OH) was found"
    
    return False, "Structure does not meet ceramide criteria"
    
# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"  # e.g. N-eicosanoylsphinganine
# result, reason = is_ceramide(test_smiles)
# print(result, reason)