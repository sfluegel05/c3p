"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
#!/usr/bin/env python
"""
Classifies: Tetradecanoate ester

A tetradecanoate ester is a fatty acid ester obtained by condensation of the carboxyl group
of tetradecanoic (myristic) acid (CH3(CH2)12CO–) with a hydroxy group of an alcohol or phenol.
To qualify, the ester group must have:
  • A carbonyl carbon (C) that is double‐bonded to an oxygen (carbonyl oxygen) and single‐bonded to:
      - an oxygen (the ester oxygen bridging to the rest of the molecule)
      - a carbon (the start of the tetradecanoate acyl chain)
  • The acyl chain (starting with the carbonyl carbon counted as carbon #1) must be linear
    (i.e. unbranched), fully saturated (all C–C bonds are single bonds) and consist exactly of 14 carbons.
Any deviation (branching, unsaturation, or chain length) disqualifies the moiety.
This function uses RDKit to search for at least one ester group meeting these criteria.
"""
from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines whether the molecule contains an ester functional group in which the acyl chain 
    is an unbranched, fully saturated chain of exactly 14 carbons (tetradecanoate).
    The method identifies candidate ester carbonyl carbons and then traverses the acyl chain,
    ensuring that every successive carbon is connected by a single bond with no branching.
    For the terminal (methyl) carbon, it checks that exactly three hydrogens are attached.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): Tuple with True and an explanation if a matching ester group is found;
                     otherwise, False and a message indicating why no matching ester was detected.
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are available.
    mol = Chem.AddHs(mol)
    
    # Iterate over all carbon atoms looking for candidate ester carbonyl carbons.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        
        carbonyl_o = None
        ester_o = None
        acyl_c = None
        
        # Examine the bonds from the candidate carbonyl carbon.
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            # Identify double bond to oxygen (carbonyl oxygen)
            if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                carbonyl_o = nbr
            # Identify single bond to oxygen (ester oxygen)
            elif nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.SINGLE:
                ester_o = nbr
            # Identify single bond to carbon (begin acyl chain)
            elif nbr.GetAtomicNum() == 6 and bond.GetBondType() == Chem.BondType.SINGLE:
                acyl_c = nbr
        
        # If we did not find the three features then skip candidate.
        if not (carbonyl_o and ester_o and acyl_c):
            continue
        
        # Now verify the acyl chain.
        # We count the carbonyl carbon as position #1.
        chain_count = 1
        prev = atom  # starting at carbonyl carbon.
        current = acyl_c
        
        valid_chain = True
        while True:
            # Verify that the bond from previous to current is a single bond (should be by selection)
            bond = mol.GetBondBetweenAtoms(prev.GetIdx(), current.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                valid_chain = False
                break

            # Ensure current atom is carbon.
            if current.GetAtomicNum() != 6:
                valid_chain = False
                break
            
            chain_count += 1
            
            # Identify next carbon(s) in the chain.
            next_carbons = []
            for nbr in current.GetNeighbors():
                # Skip hydrogens
                if nbr.GetAtomicNum() != 6:
                    continue
                # Do not go back to the previous atom.
                if nbr.GetIdx() == prev.GetIdx():
                    continue
                # Check that the connecting bond is single.
                nxt_bond = mol.GetBondBetweenAtoms(current.GetIdx(), nbr.GetIdx())
                if nxt_bond and nxt_bond.GetBondType() == Chem.BondType.SINGLE:
                    next_carbons.append(nbr)
                    
            # For a linear acyl chain, at an intermediate position we expect exactly one next carbon.
            if len(next_carbons) > 1:
                valid_chain = False  # branching detected
                break
            elif len(next_carbons) == 0:
                # We've reached the terminal carbon.
                # Confirm it is a methyl group (should have exactly three attached hydrogens).
                if current.GetTotalNumHs() != 3:
                    valid_chain = False
                break
            else:
                # Proceed to next carbon.
                prev = current
                current = next_carbons[0]
        
        # The acyl chain must have exactly 14 carbons (including the carbonyl carbon).
        if valid_chain and chain_count == 14:
            return True, ("Found an ester group with a saturated, unbranched 14‑carbon acyl chain "
                           "(tetradecanoate).")
    
    return False, "No ester group with an unbranched, saturated 14‑carbon acyl chain (tetradecanoate) was detected."


# Example usage:
# test_smiles = "CCCCCCCCCCCCCC(=O)OC"  # methyl tetradecanoate
# result, reason = is_tetradecanoate_ester(test_smiles)
# print(result, reason)