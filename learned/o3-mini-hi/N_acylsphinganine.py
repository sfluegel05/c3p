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
  3. A sphinganine backbone is present. Here we require that one of two substructure patterns is found,
     with the acyl amide nitrogen as the beginning of the sphinganine backbone:
         - Pattern without substitution on the hydroxymethyl group: "N[C;H1](CO)[C;H1](O)"
         - Pattern allowing substitution of the hydroxymethyl group (e.g. glycosylation): "N[C;H1](CO[*])[C;H1](O)"
  4. The acyl amide N in the backbone must match the one found in the C(=O)N fragment.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES.
      2. Find exactly one acyl amide bond ("C(=O)N"); if not exactly one is found, fail.
      3. Starting from the acyl carbon (the carbonyl carbon), identify a fatty acyl chain – using a breadth-first search to count contiguous carbon atoms.
      4. Look for a sphinganine backbone pattern beginning with the amide nitrogen. Two patterns are allowed:
             a. "N[C;H1](CO)[C;H1](O)" for an unsubstituted hydroxymethyl group.
             b. "N[C;H1](CO[*])[C;H1](O)" for a substituted version (e.g. glycosylated).
         In both cases the first atom (N) in the pattern should be the acyl amide nitrogen.
      5. If all tests pass, return True with an explanation.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule qualifies as an N-acylsphinganine, False otherwise.
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
    
    # Identify indices: carbonyl carbon and its bonded amide nitrogen.
    match = amide_matches[0]
    carbonyl_idx = match[0]  # carbon of the C(=O)N motif
    amide_nitrogen_idx = match[1]  # the N attached to the carbonyl
    
    # --- Step 2: Verify the fatty acyl chain from the acyl carbon ---
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    fatty_neighbor = None
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetIdx() == amide_nitrogen_idx:
            continue  # skip the nitrogen
        if neighbor.GetAtomicNum() == 6:  # look for a carbon
            fatty_neighbor = neighbor
            break
    if fatty_neighbor is None:
        return False, "Fatty acyl chain not found attached to the acyl carbon"
    
    # Use breadth-first search to count contiguous carbons in the fatty acyl chain.
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
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in chain_carbons:
                queue.append(nbr)
    if len(chain_carbons) < 8:
        return False, f"Fatty acyl chain too short; found {len(chain_carbons)} carbon(s) but expected at least 8"
    
    # --- Step 3: Look for a sphinganine backbone pattern that includes the amide nitrogen ---
    # Using less strict, non-stereospecific SMARTS.
    sphinganine_smarts_list = [
        "N[C;H1](CO)[C;H1](O)",
        "N[C;H1](CO[*])[C;H1](O)"
    ]
    sphinganine_found = False
    for smarts in sphinganine_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match_indices in matches:
            # Our pattern is defined so that the first atom is the amide nitrogen.
            if amide_nitrogen_idx not in match_indices:
                continue
            # Ensure the bond types within the backbone are single bonds.
            valid = True
            for i in range(len(match_indices) - 1):
                a_idx = match_indices[i]
                b_idx = match_indices[i + 1]
                bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
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
    
    return True, "Contains a unique acyl amide bond, a sphinganine-like (possibly glycosylated) backbone, and an adequately long fatty acyl chain"

# Example usage:
# test_smiles = "CCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"  # Example: N-(2-hydroxyeicosanoyl)sphinganine
# result, reason = is_N_acylsphinganine(test_smiles)
# print(result, reason)