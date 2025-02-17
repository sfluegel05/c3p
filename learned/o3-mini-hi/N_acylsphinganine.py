"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
#!/usr/bin/env python
"""
Classifies: N-acylsphinganine
Definition: A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
This implementation requires:
  1. Exactly one acyl amide bond (C(=O)N) in the molecule.
  2. That the acyl carbon (of the acyl amide) is attached to a fatty acyl chain (with at least 8 contiguous carbons).
  3. That a sphinganine backbone is present; here we require a match to one of our accepted substructures:
       a. The “plain” pattern: N[C@@H](CO)[C@H](O) or N[C@H](CO)[C@H](O)
       b. A variant allowing the hydroxymethyl group to be substituted (e.g. glycosylated), by using “CO[*]”
  4. That the amide nitrogen involved in the C(=O)N fragment is part of this sphinganine backbone.
Note: This approach relaxes the earlier overly strict check on the amide nitrogen’s degree.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES string.
      2. Look for exactly one acyl amide bond using the SMARTS "C(=O)N".
      3. From that bond, identify the acyl (fatty) group and verify that it contains a contiguous chain
         of at least 8 carbon atoms.
      4. Search for a sphinganine-like backbone. We allow two types of patterns:
             - One without substitution on the hydroxymethyl group:
                     "N[C@@H](CO)[C@H](O)" and "N[C@H](CO)[C@H](O)"
             - And versions where the –CH2OH group is substituted (e.g. glycosylated):
                     "N[C@@H](CO[*])[C@H](O)" and "N[C@H](CO[*])[C@H](O)"
         We also require that one of these backbone matches includes the amide nitrogen from the acyl bond.
      5. If all checks pass, return True with an explanation.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule qualifies as an N-acylsphinganine, False otherwise.
      str: Explanation detailing the decision.
    """
    # Parse the SMILES string to obtain the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Step 1: Find exactly one acyl amide bond.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} acyl amide bond(s); expected exactly 1 for N-acylsphinganine"

    # Identify the carbonyl carbon and the amide nitrogen.
    match = amide_matches[0]
    carbonyl_idx = match[0]  # carbon atom of the C(=O)N fragment
    amide_nitrogen_idx = match[1]  # nitrogen atom bonded to the carbonyl

    # (Do not enforce a strict degree check on the amide nitrogen.)
    
    # Step 2: Verify the fatty acyl chain attached to the acyl carbon.
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    fatty_neighbor = None
    for neighbor in carbonyl_atom.GetNeighbors():
        # Skip the amide nitrogen.
        if neighbor.GetIdx() == amide_nitrogen_idx:
            continue
        # Look for a carbon; assume fatty chain consists of carbons.
        if neighbor.GetAtomicNum() == 6:
            fatty_neighbor = neighbor
            break
    if fatty_neighbor is None:
        return False, "Fatty acyl chain not found attached to the acyl carbon"

    # Use a breadth-first search to count how many contiguous carbon atoms can be reached.
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

    # Step 3: Look for a sphinganine backbone pattern that includes the amide nitrogen.
    # We allow two sets of patterns: one where the hydroxymethyl group is unsubstituted and one where it is substituted.
    sphinganine_smarts_list = [
        "N[C@@H](CO)[C@H](O)",
        "N[C@H](CO)[C@H](O)",
        "N[C@@H](CO[*])[C@H](O)",
        "N[C@H](CO[*])[C@H](O)"
    ]
    sphinganine_found = False
    for smarts in sphinganine_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match_indices in matches:
            # Check that the amide nitrogen is part of this match.
            if amide_nitrogen_idx not in match_indices:
                continue
            # Ensure that bonds between matched backbone atoms are single.
            valid = True
            for i in range(len(match_indices) - 1):
                a_idx = match_indices[i]
                b_idx = match_indices[i+1]
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
        return False, "Sphinganine backbone pattern not found (or the amide N is not part of it)"
    
    # Passed all tests.
    return True, "Contains a unique acyl amide bond, a sphinganine-like (possibly glycosylated) backbone, and an adequately long fatty acyl chain"

# Example usages (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"  # N-(2-hydroxyeicosanoyl)sphinganine
# result, reason = is_N_acylsphinganine(test_smiles)
# print(result, reason)