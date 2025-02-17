"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
#!/usr/bin/env python
"""
Classifies: Tetradecanoate ester

A tetradecanoate ester is defined as a fatty acid ester obtained by condensation of the carboxy group 
of tetradecanoic acid (myristic acid; CH3(CH2)12CO–) with a hydroxy group of an alcohol or phenol.
To be accepted an ester group must have:
  • A carbonyl carbon (C) bonded to exactly three atoms:
      - one oxygen via a double bond (the carbonyl oxygen)
      - one oxygen via a single bond (the ester oxygen connecting to the alcohol/phenol side)
      - one carbon (the beginning of the fatty acyl chain)
  • The fatty acyl chain (starting with the carbonyl carbon counted as position #1) must be
    linear (no branching), fully saturated (all C–C single bonds) and have exactly 14 carbon atoms.
Any deviation in chain length, branching, or unsaturation (other than the carbonyl C=O) disqualifies the molecule.
This function uses RDKit to search for at least one ester group that meets these criteria.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines whether the molecule contains an ester functional group in which the acyl chain 
    is an unbranched, saturated chain of exactly 14 carbons (tetradecanoate).

    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): Tuple with True and an explanation if such an ester group is found;
                     otherwise, False and a message indicating why no matching ester was detected.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop through all atoms to identify potential ester carbonyl carbons.
    # The candidate carbonyl must be carbon (atomic num 6) and have exactly 3 bonds.
    # It must have:
    #   - one double bond to an oxygen (carbonyl O),
    #   - one single bond to an oxygen (ester O), and
    #   - one single bond to a carbon (acyl chain start).
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        bonds = list(atom.GetBonds())
        if len(bonds) != 3:
            continue  # not typical of a carbonyl in an ester
        carbonyl_o = None
        ester_o = None
        acyl_c = None
        
        for bond in bonds:
            nbr = bond.GetOtherAtom(atom)
            # Look for oxygen neighbors
            if nbr.GetAtomicNum() == 8:
                # Check bond type
                if bond.GetBondType() == Chem.BondType.DOUBLE and carbonyl_o is None:
                    carbonyl_o = nbr
                elif bond.GetBondType() == Chem.BondType.SINGLE and ester_o is None:
                    ester_o = nbr
            elif nbr.GetAtomicNum() == 6:
                # Only one carbon neighbor allowed (the acyl chain start)
                if acyl_c is None and bond.GetBondType() == Chem.BondType.SINGLE:
                    acyl_c = nbr
                else:
                    # More than one carbon neighbor or non-single bond: not matching expected ester fragment.
                    acyl_c = None
                    break
        
        if not (carbonyl_o and ester_o and acyl_c):
            continue  # candidate does not match the strict ester environment
        
        # Now, walk down the acyl chain.
        # We count the carbonyl carbon as position #1.
        chain_count = 1
        prev = atom   # start at the carbonyl carbon
        current = acyl_c
        
        valid_chain = True
        while True:
            # Check that the bond between prev and current is a single bond.
            bond = mol.GetBondBetweenAtoms(prev.GetIdx(), current.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                valid_chain = False
                break
            # Current atom must be carbon.
            if current.GetAtomicNum() != 6:
                valid_chain = False
                break
            
            chain_count += 1
            
            # Find the next carbon, if any.
            next_carbons = []
            for nbr in current.GetNeighbors():
                # Skip going back to the previous atom.
                if nbr.GetIdx() == prev.GetIdx():
                    continue
                # Only consider carbon neighbors connected via a single bond.
                if nbr.GetAtomicNum() == 6:
                    next_bond = mol.GetBondBetweenAtoms(current.GetIdx(), nbr.GetIdx())
                    if next_bond and next_bond.GetBondType() == Chem.BondType.SINGLE:
                        next_carbons.append(nbr)
            # To be unbranched, there can only be at most one new carbon.
            if len(next_carbons) > 1:
                valid_chain = False  # branching detected
                break
            elif len(next_carbons) == 0:
                # End of the chain reached.
                break
            else:
                # Continue walking the chain.
                prev = current
                current = next_carbons[0]
                
        # Check that the chain length is exactly 14 carbons.
        if valid_chain and chain_count == 14:
            return True, ("Found an ester group with a saturated, unbranched 14-carbon acyl chain (tetradecanoate) "
                           "derived from tetradecanoic acid.")
    return False, "No ester group with an unbranched, saturated 14-carbon acyl chain (tetradecanoate) was detected."

# Example usage:
# test_smiles = "CCCCCCCCCCCCCC(=O)OC"  # Represents methyl tetradecanoate
# result, reason = is_tetradecanoate_ester(test_smiles)
# print(result, reason)