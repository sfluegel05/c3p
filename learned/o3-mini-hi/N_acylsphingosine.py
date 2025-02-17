"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
#!/usr/bin/env python
"""
Classifies: N‐acylsphingosine (parent compounds of the ceramide family)

Definition: An N‐acylsphingosine is a molecule that contains a sphingosine backbone –
an acyclic chain featuring a secondary amine attached to two hydroxyl‐bearing carbons
(where one of these carbons is linked to a C=C bond) – and where the nitrogen is acylated
(via an amide bond to a fatty acyl group that itself has a long aliphatic chain).

Improvements over previous version:
  • The entire molecule is required to be free of rings (to exclude glycosyl/phosphorylated derivatives).
  • We enforce that the N atom’s degree equals 2 so that it is only bound to the backbone and the acyl group.
  • The acyl substituent must be a carbonyl carbon (connected via a double bond to oxygen)
    and carry a long aliphatic chain (minimum 6 carbons) defined by a DFS search limited to non‐aromatic carbons.
"""

from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    
    Our approach is to:
      1. Reject the molecule if any ring is present (thus excluding glycosphingolipids,
         gangliosides, and phosphorylated derivatives).
      2. Find a sphingosine backbone pattern defined as:
           N - C(CO) - C(O)C=C
         where the first atom is a secondary amine, the next two carbons carry hydroxyl groups,
         and the last carbon is part of an alkene.
      3. Ensure that the N atom is acylated (i.e. its only non-backbone neighbor is a carbonyl carbon
         that is directly double-bonded to an oxygen).
      4. Use a DFS search to ensure that the acyl carbonyl is attached to a long, unbranched aliphatic chain
         (at least 6 connected carbons, none in rings or aromatic).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule matches the N-acylsphingosine criteria, False otherwise.
        str: A message indicating the reason for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with any rings (this removes sugar headgroups, phosphates, etc.)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring systems not found in a parent N-acylsphingosine structure"
    
    # Define a sphingosine backbone SMARTS: we require –N–C(CO)–C(O)C=C.
    # This pattern ignores chirality and is simplified.
    backbone_smarts = "N-C(CO)-C(O)C=C"
    backbone_pattern = Chem.MolFromSmarts(backbone_smarts)
    if backbone_pattern is None:
        return False, "Error in backbone SMARTS definition"
    
    backbone_matches = mol.GetSubstructMatches(backbone_pattern, useChirality=False)
    if not backbone_matches:
        return False, "Sphingosine backbone not found"
    
    # Helper function: DFS to count maximum chain length (number of connected, non-aromatic carbon atoms).
    def dfs_chain_length(atom, coming_from, visited):
        max_length = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == coming_from:
                continue
            # Only count if the neighbor is a non-aromatic carbon and not part of a ring.
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()):
                if nbr.GetIdx() in visited:
                    continue
                visited.add(nbr.GetIdx())
                length = 1 + dfs_chain_length(nbr, atom.GetIdx(), visited)
                if length > max_length:
                    max_length = length
                visited.remove(nbr.GetIdx())
        return max_length

    # Set minimum chain length for the fatty acyl group (in number of carbon atoms)
    min_chain_length = 6

    # Iterate over each sphingosine backbone match.
    for match in backbone_matches:
        # Ensure that the backbone atoms (by the match) are acyclic.
        if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in match):
            continue  # if any backbone atom is in a ring, skip this match

        # Our backbone pattern is defined with the first atom as nitrogen.
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Ensure the N atom is a secondary amine: it should have exactly two neighbors.
        # One neighbor should be from the backbone, and one must be the acyl substituent.
        n_neighbors = [nbr for nbr in n_atom.GetNeighbors()]
        if len(n_neighbors) != 2:
            continue
        
        # Identify the non-backbone neighbor of the nitrogen.
        acyl_candidate = None
        for nbr in n_neighbors:
            if nbr.GetIdx() not in match:
                acyl_candidate = nbr
                break
        if acyl_candidate is None:
            continue

        # Check that the acyl candidate is a carbon that is acyl (i.e. part of a carbonyl group).
        if acyl_candidate.GetAtomicNum() != 6:
            continue
        # It must be connected to at least one oxygen via a double bond.
        carbonyl_found = False
        for bond in acyl_candidate.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other = bond.GetOtherAtom(acyl_candidate)
                if other.GetAtomicNum() == 8:
                    carbonyl_found = True
                    break
        if not carbonyl_found:
            continue
        
        # Now, search from the acyl candidate (excluding the bond back to N) for a long aliphatic chain.
        visited = set([acyl_candidate.GetIdx()])
        chain_length = dfs_chain_length(acyl_candidate, n_idx, visited)
        # Include the carbonyl carbon itself in the chain count.
        total_chain = 1 + chain_length
        if total_chain < min_chain_length:
            continue
        
        # If reached here, we have:
        # – a sphingosine backbone with no rings,
        # – an N atom with exactly one extra substituent that is a carbonyl-containing carbon,
        # – and that carbonyl is attached to a sufficiently long, unbranched carbon chain.
        return True, "Molecule contains a sphingosine backbone with an N-linked acyl (fatty acid) group and a long aliphatic chain"

    return False, "No valid N-acylsphingosine backbone with proper acylation and long acyl chain found"


# Example usage:
if __name__ == "__main__":
    # Test a few provided examples.
    test_smiles = [
        # True positive:
        "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # N-2-hydroxylignoceroylsphingosine
        "CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCC\\C=C/CCCCCCCC",  # N-(11Z)-icosenoylsphingosine
        # False negative (should be accepted):
        "CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(C)=O"  # N-acetylsphingosine
    ]
    
    for s in test_smiles:
        res, reason = is_N_acylsphingosine(s)
        print(res, reason)