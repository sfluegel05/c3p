"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
#!/usr/bin/env python
"""
Classifies: N-acylsphinganine
Definition: A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
This implementation requires:
  1. Exactly one acyl amide bond (C(=O)N) in the molecule.
  2. The acyl carbon is attached to a fatty acyl chain (with at least 8 contiguous carbons).
  3. A sphinganine backbone is present that begins at the acyl amide nitrogen. We allow either:
       - A pattern with an unsubstituted hydroxymethyl side chain: "N[C](CO)[C](O)"
       - Or a pattern allowing substitution (e.g. glycosylation): "N[C](CO[*])[C](O)"
  4. The acyl amide nitrogen found in the acyl bond must be part of the sphinganine backbone match.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES.
      2. Find exactly one acyl amide bond ("C(=O)N"); if not exactly one is found, classification fails.
      3. Identify the fatty acyl chain attached to the carbonyl carbon (not attached to N)
         and ensure it contains at least 8 contiguous carbon atoms.
      4. Search for a sphinganine backbone pattern that begins with the acyl amide nitrogen.
         Two patterns are allowed:
             a. "N[C](CO)[C](O)" for an unsubstituted hydroxymethyl group.
             b. "N[C](CO[*])[C](O)" for a substituted one.
         Only matches that include the acyl amide nitrogen from step 2 are accepted.
      5. If all tests pass, return True along with an explanation.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule qualifies as an N-acylsphinganine, False otherwise.
      str: Explanation detailing the decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Step 1: Find exactly one acyl amide bond ---
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} acyl amide bond(s); expected exactly 1 for N-acylsphinganine"

    # Identify the indices: carbonyl carbon and its attached amide nitrogen.
    match = amide_matches[0]
    carbonyl_idx = match[0]  # carbonyl carbon of the C(=O)N motif
    amide_nitrogen_idx = match[1]  # nitrogen in the acyl amide

    # --- Step 2: Verify the fatty acyl chain from the acyl carbon ---
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    fatty_neighbor = None
    for neighbor in carbonyl_atom.GetNeighbors():
        # Skip the neighbor that is the amide nitrogen.
        if neighbor.GetIdx() == amide_nitrogen_idx:
            continue
        # Look for a carbon attached to the carbonyl.
        if neighbor.GetAtomicNum() == 6:
            fatty_neighbor = neighbor
            break
    if fatty_neighbor is None:
        return False, "Fatty acyl chain not found attached to the acyl carbon"

    # Use breadth-first search (BFS) to count contiguous carbon atoms in the fatty acyl chain.
    visited = set()
    queue = deque([fatty_neighbor])
    chain_carbons = set()
    while queue:
        atom = queue.popleft()
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetIdx() in chain_carbons:
            continue
        chain_carbons.add(atom.GetIdx())
        # Traverse neighbors that are carbons.
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in chain_carbons:
                queue.append(nbr)
    if len(chain_carbons) < 8:
        return False, f"Fatty acyl chain too short; found {len(chain_carbons)} carbon(s) but expected at least 8"

    # --- Step 3: Look for a sphinganine backbone pattern that includes the acyl amide nitrogen ---
    # Relaxed SMARTS patterns (removing explicit hydrogen count and chirality) to match sphinganine backbone.
    sphinganine_patterns = [
        "N[C](CO)[C](O)",       # for unsubstituted hydroxymethyl side chain
        "N[C](CO[*])[C](O)"      # for substituted (e.g. glycosylated) version
    ]
    sphinganine_found = False
    for smarts in sphinganine_patterns:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match_indices in matches:
            # Ensure the acyl amide nitrogen is part of the backbone pattern match.
            if amide_nitrogen_idx not in match_indices:
                continue
            # Optionally, we could check that the bonds between the matched atoms are single bonds.
            valid = True
            for i in range(len(match_indices)-1):
                bond = mol.GetBondBetweenAtoms(match_indices[i], match_indices[i+1])
                if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                    valid = False
                    break
            if valid:
                sphinganine_found = True
                break
        if sphinganine_found:
            break
    if not sphinganine_found:
        return False, "Sphinganine backbone pattern not found (or the acyl amide N is not part of it)"

    return True, "Contains a unique acyl amide bond, an appropriate sphinganine backbone, and an adequately long fatty acyl chain"