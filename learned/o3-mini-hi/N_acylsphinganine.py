"""
Classifies: CHEBI:31488 N-acylsphinganine
"""
#!/usr/bin/env python
"""
Classifies: N-acylsphinganine
Definition: A ceramide consisting of sphinganine in which one of the amino hydrogens is substituted by a fatty acyl group.
This implementation requires:
  1. Exactly one acyl amide bond (C(=O)N) in the molecule.
  2. The acyl carbon must be attached to a fatty acyl chain of at least 8 contiguous carbons.
  3. A sphinganine backbone that begins at the acyl amide nitrogen is detected. Two SMARTS patterns are used:
       a. "N[C](CO)[C](O)" for the unsubstituted hydroxymethyl version.
       b. "N[C](CO[*])[C](O)" to allow for substitution (e.g. glycosylation).
     In either case the match must include the acyl amide N.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_N_acylsphinganine(smiles: str):
    """
    Determines if a molecule is an N-acylsphinganine based on its SMILES string.
    
    Strategy:
      1. Parse the SMILES.
      2. Find exactly one acyl amide bond ("C(=O)N").
      3. Identify the fatty acyl chain attached to the acyl carbon (the one not attached to N)
         and ensure it contains at least 8 contiguous carbon atoms.
      4. Search for a sphinganine backbone starting at the acyl amide nitrogen.
         We use two SMARTS patterns on a copy of the molecule from which stereochemistry is removed,
         and require that the acyl amide nitrogen is part of the match.
      5. If all tests pass, return True with an explanation.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      bool: True if the molecule qualifies as an N-acylsphinganine, False otherwise.
      str: Explanation detailing the decision.
    """
    # --- Step 1: Parse the SMILES string ---
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # --- Step 2: Find exactly one acyl amide bond (C(=O)N) ---
    amide_pattern = Chem.MolFromSmarts("C(=O)N")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) != 1:
        return False, f"Found {len(amide_matches)} acyl amide bond(s); expected exactly 1 for N-acylsphinganine"

    # Record the indices for the acyl amide match:
    # match[0] = carbonyl carbon, match[1] = amide nitrogen.
    match = amide_matches[0]
    carbonyl_idx = match[0]
    amide_nitrogen_idx = match[1]

    # --- Step 3: Verify the fatty acyl chain attached to the carbonyl carbon ---
    carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
    fatty_neighbor = None
    # Find a neighbor of the carbonyl that is not the amide nitrogen.
    for nbr in carbonyl_atom.GetNeighbors():
        if nbr.GetIdx() == amide_nitrogen_idx:
            continue
        if nbr.GetAtomicNum() == 6:  # carbon
            fatty_neighbor = nbr
            break
    if fatty_neighbor is None:
        return False, "Fatty acyl chain not found attached to the acyl carbon"

    # Now count contiguous carbon atoms starting from this neighbor.
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
        return False, f"Fatty acyl chain too short; found {len(chain_carbons)} carbon(s) but at least 8 are required"

    # --- Step 4: Look for a sphinganine backbone pattern that includes the acyl amide nitrogen ---
    # To avoid stereochemistry issues, work on a copy of the molecule with chirality information removed.
    mol_no_stereo = Chem.Mol(mol)
    Chem.RemoveStereochemistry(mol_no_stereo)
    
    # Two patterns are allowed to match the sphinganine backbone.
    sphinganine_smarts = [
        "N[C](CO)[C](O)",       # unsubstituted hydroxymethyl side chain
        "N[C](CO[*])[C](O)"      # substituted (e.g. glycosylated) version
    ]
    backbone_found = False
    for smart in sphinganine_smarts:
        pattern = Chem.MolFromSmarts(smart)
        if pattern is None:
            continue
        matches = mol_no_stereo.GetSubstructMatches(pattern)
        for match_indices in matches:
            # Ensure that the acyl amide nitrogen is part of this backbone fragment.
            if amide_nitrogen_idx in match_indices:
                # Optionally check that bonds in the found pattern are single bonds.
                valid = True
                for i in range(len(match_indices)-1):
                    bond = mol_no_stereo.GetBondBetweenAtoms(match_indices[i], match_indices[i+1])
                    if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                        valid = False
                        break
                if valid:
                    backbone_found = True
                    break
        if backbone_found:
            break

    if not backbone_found:
        return False, "Sphinganine backbone pattern not found (or the acyl amide N is not part of it)"
    
    return True, "Contains a unique acyl amide bond, an adequately long fatty acyl chain, and a sphinganine backbone starting at the acyl amide nitrogen"

# (Optional) For simple testing uncomment the following lines:
# test_smiles = "CCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)CCCCCCCCCCCCCCC"
# result, reason = is_N_acylsphinganine(test_smiles)
# print(result, reason)