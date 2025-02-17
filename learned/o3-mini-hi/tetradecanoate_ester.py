"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
#!/usr/bin/env python
"""
Classifies: Tetradecanoate ester

A tetradecanoate ester is defined as a fatty acid ester obtained by condensation of the
carboxyl group of tetradecanoic (myristic) acid (CH3(CH2)12CO–) with a hydroxy group of
an alcohol or phenol. To be accepted an ester group must have:
  • A carbonyl carbon (C) bonded to exactly three heavy atoms:
      - one oxygen via a double bond (the carbonyl oxygen)
      - one oxygen via a single bond (the ester oxygen linking to the rest of the molecule)
      - one carbon (the beginning of the fatty acyl chain)
  • The acyl chain (starting with the carbonyl carbon counted as carbon #1) must be linear 
    (i.e. unbranched), fully saturated (only C–C single bonds), and have exactly 14 carbons.
Any deviation in chain length, branching, or unsaturation (aside from the carbonyl C=O) disqualifies the molecule.
This function uses RDKit to search for at least one ester group that meets these tight criteria.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines whether the molecule contains an ester functional group in which the acyl chain 
    is an unbranched, saturated chain of exactly 14 carbons (tetradecanoate).
    It ensures that the candidate chain is truly linear by adding hydrogens and by ensuring that the terminal 
    carbon is a methyl group (which should have exactly three hydrogens).

    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): Tuple with True and an explanation if such an ester group is found;
                     otherwise, False and a message indicating why no matching ester was detected.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so that hydrogen counts are available.
    mol = Chem.AddHs(mol)
    
    # Loop through all atoms to identify candidate ester carbonyl carbons.
    # We look only at heavy atom neighbors (non-hydrogen) for degree counts.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Only consider heavy atoms (non-H) when checking the bonding environment.
        heavy_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 3:
            continue
        
        # For a correct ester carbonyl carbon we need:
        #   - exactly one double bond to an oxygen (carbonyl oxygen)
        #   - exactly one single bond to an oxygen (ester oxygen)
        #   - exactly one single bond to a carbon (start of the acyl chain)
        carbonyl_o = None
        ester_o = None
        acyl_c = None
        
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.BondType.DOUBLE and carbonyl_o is None:
                    carbonyl_o = nbr
                elif bond.GetBondType() == Chem.BondType.SINGLE and ester_o is None:
                    ester_o = nbr
            elif nbr.GetAtomicNum() == 6:
                # Only accept a single carbon neighbor (for the acyl chain) via a single bond.
                if bond.GetBondType() == Chem.BondType.SINGLE and acyl_c is None:
                    acyl_c = nbr
                else:
                    # If there is more than one carbon neighbor or non-single bond, skip candidate.
                    acyl_c = None
                    break
        
        if not (carbonyl_o and ester_o and acyl_c):
            continue  # candidate does not match strict ester environment
        
        # Now, traverse the acyl chain.
        # Count the carbonyl carbon as position #1.
        chain_count = 1
        prev = atom   # start at the carbonyl carbon (position 1)
        current = acyl_c
        
        valid_chain = True
        while True:
            # Ensure bond between prev and current is a single bond.
            bond = mol.GetBondBetweenAtoms(prev.GetIdx(), current.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                valid_chain = False
                break
            
            # The current atom must be a carbon.
            if current.GetAtomicNum() != 6:
                valid_chain = False
                break
            
            chain_count += 1
            
            # Look for next carbon in the chain.
            next_carbons = []
            for nbr in current.GetNeighbors():
                # Limit to heavy atoms only.
                if nbr.GetAtomicNum() == 1:
                    continue
                # Avoid going back to the previous atom.
                if nbr.GetIdx() == prev.GetIdx():
                    continue
                # Must be connected by a single bond.
                nxt_bond = mol.GetBondBetweenAtoms(current.GetIdx(), nbr.GetIdx())
                if nxt_bond and nxt_bond.GetBondType() == Chem.BondType.SINGLE and nbr.GetAtomicNum() == 6:
                    next_carbons.append(nbr)
            # For an unbranched chain, there should be at most one new carbon.
            if len(next_carbons) > 1:
                valid_chain = False  # branching detected
                break
            elif len(next_carbons) == 0:
                # End of chain; now check that the terminal carbon is a methyl i.e. has 3 hydrogens.
                if current.GetTotalNumHs() != 3:
                    valid_chain = False
                break
            else:
                prev = current
                current = next_carbons[0]
                
        # The acyl chain must be exactly 14 carbons long including the carbonyl carbon.
        if valid_chain and chain_count == 14:
            return True, ("Found an ester group with a saturated, unbranched 14‑carbon acyl chain "
                           "(tetradecanoate) derived from tetradecanoic acid.")
    return False, "No ester group with an unbranched, saturated 14‑carbon acyl chain (tetradecanoate) was detected."


# Example usage:
# test_smiles = "CCCCCCCCCCCCCC(=O)OC"  # methyl tetradecanoate
# result, reason = is_tetradecanoate_ester(test_smiles)
# print(result, reason)