"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
#!/usr/bin/env python
"""
Classifies: N-acylsphinganine 
Definition: A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
This function attempts to identify such molecules by requiring a single acyl amide bond (C(=O)N)
and a sphinganine backbone (identified roughly by the substructure N[C@@H](CO)[C@H](O) or its alternative stereochemistry).
It also performs a simple check on the length of the fatty acyl chain.
"""

from rdkit import Chem
from collections import deque

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is N-acylsphinganine based on its SMILES string.
    
    The strategy is as follows:
     1. Parse the SMILES.
     2. Ensure there is exactly one amide bond (identified by the substructure C(=O)N).
     3. Check that the region adjacent to the nitrogen atom contains a sphinganine-like backbone,
        as identified by common substructure patterns (with stereochemical markers).
     4. Do a simple check that the acyl (fatty) chain attached via the carbonyl is long enough 
        (here we require at least 8 contiguous carbon atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if molecule is N-acylsphinganine, False otherwise.
        str: Reason for the classification.
    """
    # Step 1: Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Step 2: Look for exactly one acyl amide bond.
    # The SMARTS "C(=O)N" identifies a carbonyl carbon directly attached to a nitrogen.
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} amide bond(s); expected exactly 1 for N-acylsphinganine"
    
    # Step 3: Check for the sphinganine backbone.
    # Many sphinganine derivatives include a backbone fragment like: N[C@@H](CO)[C@H](O) 
    # (or with alternate stereochemistry). We will try both.
    sphinganine_smarts_list = [
        "N[C@@H](CO)[C@H](O)",
        "N[C@H](CO)[C@H](O)"
    ]
    sphinganine_found = False
    for smarts in sphinganine_smarts_list:
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        if mol.HasSubstructMatch(pattern):
            sphinganine_found = True
            break
    if not sphinganine_found:
        return False, "Sphinganine backbone pattern not found"
    
    # Step 4: Check that the acyl (fatty) chain is sufficiently long.
    # In an N-acylsphinganine the fatty acyl group is attached via the amide bond.
    # We already got the amide match – the pattern "C(=O)N" returns a tuple (carbonyl_carbon, 
    # amide_nitrogen). The fatty acyl chain is attached on the carbonyl side (the acyl carbon),
    # excluding the neighbor that came from the N.
    match = amide_matches[0]
    carbonyl_idx = match[0]
    amide_nitrogen_idx = match[1]
    
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    fatty_neighbor = None
    for neighbor in carbonyl_atom.GetNeighbors():
        if neighbor.GetIdx() == amide_nitrogen_idx:
            continue
        if neighbor.GetAtomicNum() == 6:  # looking for a carbon
            fatty_neighbor = neighbor
            break
    if fatty_neighbor is None:
        return False, "Fatty acyl chain not found attached to the acyl carbon"
    
    # We now perform a breadth-first search (restricted to carbon atoms) starting from fatty_neighbor
    # to count the number of contiguous carbons in this chain.
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
    
    # Passed all tests, we classify as an N-acylsphinganine.
    return True, "Contains a single acyl amide bond with a sphinganine backbone and adequate fatty acyl chain length"

# Example usage (uncomment to test):
# test_smiles = "CCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"
# result, reason = is_N_acylsphinganine(test_smiles)
# print(result, reason)