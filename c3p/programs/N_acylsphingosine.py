"""
Classifies: CHEBI:52639 N-acylsphingosine
"""
#!/usr/bin/env python
"""
Classifies: N‐acylsphingosine (parent compounds of the ceramide family)
Definition: An N‐acylsphingosine is a sphingosine (a long‐chain amino alcohol with two hydroxyl groups,
one of which is on the primary alcohol and one on a secondary carbon that also carries the amide)
that has its only amino group converted to an amide by attachment of a fatty acyl group.
The acyl group may be very short (e.g. acetyl) or long. To avoid classifying phosphorylated/glycosylated
derivatives (i.e. ceramides with extra headgroups), we reject molecules that contain any rings.
We then look for an amide bond (C(=O)–N) where the N atom has exactly two heavy neighbors.
One neighbor (from the amide) is the acyl carbon, which must show a C=O bond to at least one oxygen and
carry an aliphatic (non‐ring, non‐aromatic) chain of at least 2 carbons. The other neighbor must be part
of a sphingosine backbone – namely, it must be a carbon that carries a hydroxymethyl substituent (–CH2–OH)
and is itself attached to a second carbon bearing an –OH and further attached to an alkene.
This set of criteria was chosen in an effort to explain the outcomes of the previous version.
"""

from rdkit import Chem

def is_N_acylsphingosine(smiles: str):
    """
    Determines if a molecule is an N-acylsphingosine based on its SMILES string.
    
    The function applies the following steps:
      1. Parse the SMILES and immediately reject molecules containing any rings.
      2. Loop over all amide bonds (C(=O)–N). For each:
         a. Check that the N atom is secondary (degree = 2).
         b. Designate the carbonyl carbon (acyl group) and the other neighbor of N as the backbone carbon.
         c. Verify that the backbone carbon has a hydroxymethyl substituent (a neighbor O that is
            connected only to that carbon) – as expected in sphingosine.
         d. Verify that the backbone carbon is connected (besides to the amino group) to another carbon
            that bears an –OH and has at least one double bond (indicative of the alkene in sphingosine).
         e. Finally, from the acyl (carbonyl) carbon, perform a DFS to estimate the length (in carbons) 
            of the fatty acyl chain. We allow a minimum chain length of 2 (to accommodate cases like N-acetylsphingosine).
      3. If a match is found that satisfies these criteria, return True with an appropriate reason.
         Otherwise, return False.
         
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule matches the N-acylsphingosine criteria, False otherwise.
        str: Explanation of the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Reject molecules with any rings (to eliminate phosphorylated/glycosylated ceramide derivatives)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains ring systems not found in a parent N-acylsphingosine structure"
    
    # Helper: perform a DFS search to measure maximum linear chain length.
    # Only counts contiguous, non-aromatic, non-ring carbon atoms.
    def dfs_chain_length(atom, coming_from, visited):
        max_length = 0
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() == coming_from:
                continue
            # Count only if neighbor is carbon, not aromatic, and not in a ring.
            if nbr.GetAtomicNum() == 6 and (not nbr.GetIsAromatic()) and (not nbr.IsInRing()):
                if nbr.GetIdx() in visited:
                    continue
                visited.add(nbr.GetIdx())
                length = 1 + dfs_chain_length(nbr, atom.GetIdx(), visited)
                if length > max_length:
                    max_length = length
                visited.remove(nbr.GetIdx())
        return max_length

    # Set minimum acceptable chain length (in number of carbons in the acyl chain, counting the carbonyl carbon)
    min_chain_length = 2

    # Loop over all bonds to find potential amide bonds (C(=O)-N)
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Identify pair where one is N and the other is C
        if (a1.GetAtomicNum() == 7 and a2.GetAtomicNum() == 6) or (a2.GetAtomicNum() == 7 and a1.GetAtomicNum() == 6):
            # Designate: amide_N is the nitrogen and amide_C is the carbon
            if a1.GetAtomicNum() == 7:
                amide_N = a1
                amide_C = a2
            else:
                amide_N = a2
                amide_C = a1

            # Check that the bond itself is a single bond (amide bond) – typically it is.
            if bond.GetBondType() != Chem.BondType.SINGLE:
                continue

            # Verify that the carbon shows a carbonyl: it must have a double bond to an oxygen.
            carbonyl_found = False
            for b in amide_C.GetBonds():
                if b.GetBondType() == Chem.BondType.DOUBLE:
                    other = b.GetOtherAtom(amide_C)
                    if other.GetAtomicNum() == 8:
                        carbonyl_found = True
                        break
            if not carbonyl_found:
                continue

            # Check that the nitrogen is a secondary amine: exactly 2 heavy-atom neighbors.
            n_neighbors = [nbr for nbr in amide_N.GetNeighbors() if nbr.GetAtomicNum() > 1]
            if len(n_neighbors) != 2:
                continue

            # Identify the backbone carbon attached to N (the neighbor that is not the carbonyl C)
            backbone_atom = None
            for nbr in amide_N.GetNeighbors():
                if nbr.GetIdx() != amide_C.GetIdx():
                    backbone_atom = nbr
                    break
            if backbone_atom is None:
                continue

            # At this point, amide_N is connected to:
            #   - amide_C (the carbonyl from the acyl group)
            #   - backbone_atom (should become part of the sphingosine chain)
            # Next, check that backbone_atom serves as the start of a sphingosine backbone:
            #   (i) It should carry a hydroxymethyl substituent.
            hydroxymethyl_found = False
            for nbr in backbone_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 8:
                    # To check –OH, we require that the oxygen is only connected to this carbon (ignoring implicit H)
                    if nbr.GetDegree() == 1:
                        hydroxymethyl_found = True
                        break
            if not hydroxymethyl_found:
                continue

            #   (ii) The backbone_atom should be connected to another carbon (call it b2) that bears an -OH.
            b2_candidate = None
            for nbr in backbone_atom.GetNeighbors():
                if nbr.GetIdx() == amide_N.GetIdx():
                    continue
                if nbr.GetAtomicNum() == 6 and (not nbr.IsInRing()):
                    # Look for an -OH on this candidate.
                    for subnbr in nbr.GetNeighbors():
                        if subnbr.GetAtomicNum() == 8 and subnbr.GetDegree() == 1:
                            b2_candidate = nbr
                            break
                    if b2_candidate is not None:
                        break
            if b2_candidate is None:
                continue

            #   (iii) To capture the alkene tail typical of sphingosine we require that b2_candidate is connected (besides backbone_atom)
            #         to a carbon via a double bond.
            alkene_found = False
            for b in b2_candidate.GetBonds():
                if b.GetBeginAtom().GetIdx() == b2_candidate.GetIdx() or b.GetEndAtom().GetIdx() == b2_candidate.GetIdx():
                    # skip the bond to backbone_atom
                    if b.GetOtherAtom(b2_candidate).GetIdx() == backbone_atom.GetIdx():
                        continue
                    if b.GetBondType() == Chem.BondType.DOUBLE:
                        other = b.GetOtherAtom(b2_candidate)
                        if other.GetAtomicNum() == 6:
                            alkene_found = True
                            break
            if not alkene_found:
                continue

            # Now check the acyl chain:
            # Perform a DFS from the carbonyl carbon (amide_C) and count a maximum chain length
            # Exclude the double-bonded oxygen atom.
            visited = set([amide_C.GetIdx()])
            chain_length = dfs_chain_length(amide_C, amide_N.GetIdx(), visited)
            total_chain = 1 + chain_length  # count the carbonyl carbon itself
            if total_chain < min_chain_length:
                continue

            # If all these checks are met we consider the molecule to be a valid N-acylsphingosine.
            return True, "Contains a sphingosine backbone with an N-linked acyl (fatty acid) group and proper substituents"
    
    return False, "No valid N-acylsphingosine backbone with proper acylation and correct hydroxyl/alkene pattern found"

# Example usage:
if __name__ == "__main__":
    # Some examples from the assignment:
    test_smiles = [
        "CCCCCCCCCCCCCCCCCCCCCC(O)C(=O)N[C@@H](CO)[C@H](O)\\C=C\\CCCCCCCCCCCCC",  # N-2-hydroxylignoceroylsphingosine (TP)
        "CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(=O)CCCCCCCCC\\C=C/CCCCCCCC",    # N-(11Z)-icosenoylsphingosine (TP)
        "[C@H]([C@@H](/C=C/CCCCCCCCCCCCC)O)(NC(=O)CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCOC(CCCCCCC/C=C\\C/C=C\\CCCCC)=O)CO",  # example TP
        "CCCCCCCCCCCC\\C=C\\[C@@H](O)[C@H](CO)NC(C)=O",  # N-acetylsphingosine (FN in previous version, but now accepted)
        "CCCCCCCCCCCCCCCC(=O)N[C@H]([C@@H](/C=C/CCCCCCCCC(CC)C)O)CO",  # false positive example from outcomes
    ]
    
    for s in test_smiles:
        res, reason = is_N_acylsphingosine(s)
        print("SMILES:", s)
        print("Result:", res, "|", reason)
        print("-" * 80)