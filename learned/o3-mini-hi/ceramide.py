"""
Classifies: CHEBI:17761 ceramide
"""
"""
Classifies: Ceramide (N-acyl-sphingoid bases)

Ceramides are sphingoid base derivatives that contain an amide‐linked fatty acid.
The fatty acid is typically saturated or monounsaturated with chain lengths between 14 and 26 carbons.
Many sphingoid bases have an amino group (acylated to give the amide) and show at least one –OH (often a CH2OH group) attached to a carbon.
This implementation improves on the previous attempt by:
  – Checking for an amide group in which the carbonyl-side has a contiguous non‐cyclic alkyl chain (14–26 carbons)
  – Requiring that the amide nitrogen (base‐side) is bonded to at least one carbon that carries a hydroxyl group.
Note: This is still heuristic and may mis‐classify borderline or unusually decorated structures.
"""
from rdkit import Chem

def is_ceramide(smiles: str):
    """
    Determines if a molecule is a ceramide, based on its SMILES string.

    It:
     1. Finds at least one amide group (pattern "C(=O)N").
     2. Inspects each amide group to see whether on the carbonyl (acyl) side
        there is a contiguous (acyclic) carbon chain of length 14 to 26.
     3. Checks that on the nitrogen (base) side, the N is connected to a carbon
        that has at least one directly attached hydroxyl (–OH) group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is classified as a ceramide, False otherwise.
        str: Explanation of the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # --- 1. Look for an amide group ---
    # A simple SMARTS for an amide; here match returns (carbonyl carbon, amide nitrogen)
    amide_smarts = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_smarts)
    if not amide_matches:
        return False, "No amide group (C(=O)N) found"
    
    # Helper function to measure a contiguous, acyclic carbon chain length.
    def linear_chain_length(mol, current_idx, coming_from_idx, visited):
        """
        Recursively determine the length of the contiguous (acyclic, non-ring) carbon chain.
        Only follows atoms with atomic number 6 (carbon) and which are not in any ring.
        This DFS returns the maximum chain length (counting the starting atom as length 1).
        """
        current_atom = mol.GetAtomWithIdx(current_idx)
        # Do not continue if the atom is in a ring.
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

    # Helper function to test whether the amide N is attached to a potential sphingoid base fragment.
    def has_sphingoid_feature(n_atom, carbonyl_idx):
        """
        Check that the amide nitrogen (n_atom) is connected (aside from the carbonyl carbon)
        to at least one carbon that in turn directly bears a hydroxyl group.
        """
        for nb in n_atom.GetNeighbors():
            # Skip the carbonyl carbon
            if nb.GetIdx() == carbonyl_idx:
                continue
            if nb.GetAtomicNum() == 6:  # carbon neighbor
                # Check if this neighbor has an -OH group.
                for sub_nb in nb.GetNeighbors():
                    # Look for an attached oxygen (non-carbonyl oxygen)
                    if sub_nb.GetAtomicNum() == 8:
                        # Exclude carbonyl oxygen (which would be attached to the carbonyl carbon)
                        if sub_nb.GetDegree() == 1:
                            return True
        return False

    fatty_acyl_found = False
    sphingo_found = False
    reasons = []
    
    # --- 2 & 3. For each amide match, check the fatty acyl chain on the acyl side and sphingoid feature on the base side.
    for match in amide_matches:
        carbonyl_idx = match[0]
        n_idx = match[1]
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # For the acyl side, examine neighbors of the carbonyl carbon other than the amide nitrogen
        candidate_chain_idxs = []
        for nb in carbonyl_atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            # skip the nitrogen and the typical carbonyl oxygen
            if nb_idx == n_idx:
                continue
            if nb.GetAtomicNum() == 8:
                continue
            if nb.GetAtomicNum() == 6:
                candidate_chain_idxs.append(nb_idx)
        # Check each candidate chain
        for cand in candidate_chain_idxs:
            # Only consider candidate carbons that are not in a ring
            if mol.GetAtomWithIdx(cand).IsInRing():
                continue
            length = linear_chain_length(mol, cand, carbonyl_idx, {cand})
            # Check if the chain length is within typical fatty acid range (14 – 26 carbons).
            if 14 <= length <= 26:
                fatty_acyl_found = True
                break  # We have found a valid fatty acyl chain for this amide.
        if not fatty_acyl_found:
            continue  # Try the next amide match
        
        # Next, check the sphingoid base: the amide nitrogen should have a neighbor that is a carbon
        # with a directly attached hydroxyl group.
        if has_sphingoid_feature(n_atom, carbonyl_idx):
            sphingo_found = True
            # If both required features are present, we can return True.
            return True, ("Molecule contains an amide-linked fatty acid (chain length ~14–26 carbons) "
                          "and a sphingoid base feature (N bonded to a carbon bearing an –OH group)")
        else:
            reasons.append("Amide group found with fatty acyl chain but no sphingoid base feature "
                           "(N not bonded to a carbon with –OH) for this linkage.")
    
    if not fatty_acyl_found:
        return False, "No fatty acyl chain with 14–26 contiguous carbons (attached to an amide carbonyl) found"
    if not sphingo_found:
        # If we did find an amide with fatty acyl but none of them passed the sphingoid base test
        return False, "Fatty acyl chain detected, but no suitable sphingoid base feature (N-adjacent –OH) was found"
    
    return False, "Unknown structure"
    
# Example usage (uncomment to test):
# smiles_test = "CCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"  # e.g. N-eicosanoylsphinganine
# result, reason = is_ceramide(smiles_test)
# print(result, reason)