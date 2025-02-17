"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
#!/usr/bin/env python
"""
Classifies: N-acylsphinganine 
Definition: A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
This function identifies such molecules by:
  1. Ensuring exactly one acyl amide bond (C(=O)N).
  2. Checking that the acyl (fatty) chain attached via the carbonyl is long (at least 8 contiguous carbons).
  3. Looking for a sphinganine backbone fragment. Here we require a match to patterns such as
     "N[C@@H](CO)[C@H](O)" (or an alternative stereochemistry) that is fully saturated (i.e. no C=C bonds along the backbone).
  4. Verifying that the amide nitrogen is only connected to the carbonyl carbon and the sphinganine backbone.
Note: Some false positives previously arose because unsaturated sphingoid backbones (sphingosines) or extra substituents were present.
"""

from rdkit import Chem
from collections import deque

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES.
      2. Ensure there is exactly one acyl amide bond (pattern "C(=O)N").
      3. Check that the amide nitrogen has exactly two heavy-atom neighbors (one from the acyl group and one from the sphinganine backbone).
      4. Look for a sphinganine backbone fragment, requiring a match to accepted patterns and that the bonds between atoms in the substructure are single.
      5. Perform a breadth-first search (restricted to carbon atoms) from the acyl carbon (on the fatty side) to count the number of contiguous carbons.
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if molecule qualifies as an N-acylsphinganine, False otherwise.
      str: Explanation detailing the decision.
    """
    # Step 1: Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2: Search for exactly one acyl amide bond. 
    # The SMARTS "C(=O)N" identifies a carbonyl carbon directly attached to a nitrogen.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide bond(s); expected exactly 1 for N-acylsphinganine"
    
    # Get the specific match and identify the indices.
    match = amide_matches[0]
    carbonyl_idx = match[0]  # carbon of C(=O)N
    amide_nitrogen_idx = match[1]  # nitrogen attached to the carbonyl
    
    # Extra check on the amide nitrogen connectivity:
    amide_nitrogen = mol.GetAtomWithIdx(amide_nitrogen_idx)
    # In a valid N-acylsphinganine the nitrogen should be originally NH2, and after acylation should have:
    #  - one bond to the carbonyl carbon and one bond to the sphinganine backbone.
    # (Hydrogens are implicit so heavy-atom degree should be exactly 2.)
    if amide_nitrogen.GetDegree() != 2:
        return False, "Amide nitrogen has extra substituents; not a valid N-acylsphinganine"
    
    # Step 3: Look for a sphinganine backbone.
    # We test a couple of stereochemical variants.
    sphinganine_smarts_list = [
        "N[C@@H](CO)[C@H](O)",
        "N[C@H](CO)[C@H](O)"
    ]
    sphinganine_found = False
    for smarts in sphinganine_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match_indices in matches:
            # For the match to represent a sphinganine backbone, we require that 
            # all bonds between the backbone atoms are single (no double bonds).
            valid = True
            # Get the atoms in the match
            for i in range(len(match_indices)-1):
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
        return False, "Sphinganine (saturated backbone) pattern not found"
    
    # Step 4: Check that the fatty acyl chain is sufficiently long.
    # The acyl chain is connected to the carbonyl carbon (other than via the amide bond).
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    fatty_neighbor = None
    for neighbor in carbonyl_atom.GetNeighbors():
        # Skip the amide nitrogen
        if neighbor.GetIdx() == amide_nitrogen_idx:
            continue
        if neighbor.GetAtomicNum() == 6:  # looking for a carbon atom
            fatty_neighbor = neighbor
            break
    if fatty_neighbor is None:
        return False, "Fatty acyl chain not found attached to the acyl carbon"
    
    # Count contiguous carbon atoms (allowing branches) in the fatty chain using a breadth-first search.
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
    
    # Passed all tests.
    return True, "Contains a single acyl amide bond with a saturated sphinganine backbone and an adequately long fatty acyl chain"

# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"  # N-(2-hydroxyeicosanoyl)sphinganine
# result, reason = is_N_acylsphinganine(test_smiles)
# print(result, reason)