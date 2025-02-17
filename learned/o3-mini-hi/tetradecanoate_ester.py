"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
#!/usr/bin/env python
"""
Classifies: Tetradecanoate ester

A tetradecanoate ester is defined as a fatty acid ester obtained by condensation of the carboxyl group 
of tetradecanoic acid (myristic acid; CH3(CH2)12CO–) with a hydroxy group of an alcohol or phenol.
Thus, to be accepted an ester group must have:
  • A carbonyl carbon (C) bonded to exactly:
      - one oxygen via a double bond (the carbonyl oxygen)
      - one oxygen via a single bond (the ester oxygen that connects to the alcohol/phenol residue)
      - one carbon (the start of the acyl chain)
  • A linear, unbranched acyl chain (traversed via only single bonds and only through carbon atoms)
    starting at the acyl neighbor such that when “walking” the chain (counting the carbonyl as 1)
    exactly 14 carbons are encountered.
  • Any unsaturation (other than the carbonyl C=O), branching, or acyl chain length not equal to 14
    will rule out a match.
This code uses RDKit to look for at least one ester group that meets these criteria.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines whether the molecule contains an ester functional group in which the acyl chain 
    is an unbranched, saturated chain of exactly 14 carbons (tetradecanoate), derived from tetradecanoic acid.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): Tuple with True and a matching reason if such an ester group is found;
                     otherwise, False and a message indicating why no such group was detected.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop over all atoms and check for candidate ester carbonyl carbons
    for atom in mol.GetAtoms():
        # Only consider carbon atoms
        if atom.GetAtomicNum() != 6:
            continue

        # Look at bonds to see if there is one double-bond to an oxygen (carbonyl oxygen)
        # and one single bond to an oxygen (ester oxygen), plus exactly one carbon neighbor.
        carbonyl_oxygen = None
        ester_oxygen = None
        acyl_neighbor = None   # the carbon that begins the fatty acyl chain
        
        bonds = atom.GetBonds()
        # We also expect that the candidate carbonyl forms part of a –C(=O)–O– fragment.
        for bond in bonds:
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                # Look at the bond order
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    if carbonyl_oxygen is None:
                        carbonyl_oxygen = nbr
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    if ester_oxygen is None:
                        ester_oxygen = nbr
            elif nbr.GetAtomicNum() == 6:
                # This should be the acyl chain start; make sure there is exactly one.
                if acyl_neighbor is None:
                    acyl_neighbor = nbr
                else:
                    # More than one carbon neighbor gives ambiguity
                    acyl_neighbor = None
                    break  # break out of bond loop
        
        # Check that we have one carbonyl oxygen, one ester oxygen, and exactly one acyl-chain carbon.
        if not (carbonyl_oxygen and ester_oxygen and acyl_neighbor):
            continue  # move to next candidate
        
        # Now traverse the acyl chain.
        # We count the carbonyl carbon as #1.
        chain_count = 1
        prev_atom = atom  # start with the carbonyl carbon
        current_atom = acyl_neighbor
        
        valid_chain = True
        
        # We now “walk” along the chain via single, C–C bonds.
        while True:
            # Verify that the bond from prev_atom to current_atom is a single bond.
            bond = mol.GetBondBetweenAtoms(prev_atom.GetIdx(), current_atom.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                valid_chain = False
                break
            
            # Check that current_atom is carbon (it must be, because we are looking at the acyl chain)
            if current_atom.GetAtomicNum() != 6:
                valid_chain = False
                break
            
            # For all bonds attached to current_atom (other than back to prev_atom), 
            # make sure they are single bonds (to rule out unsaturation in the chain).
            for nbr_bond in current_atom.GetBonds():
                if nbr_bond.GetBeginAtomIdx() == prev_atom.GetIdx() or nbr_bond.GetEndAtomIdx() == prev_atom.GetIdx():
                    continue
                if nbr_bond.GetBondType() != Chem.BondType.SINGLE:
                    valid_chain = False
                    break
            if not valid_chain:
                break
            
            chain_count += 1
            
            # Look for the "next" carbon in the chain.
            next_carbons = []
            for nbr in current_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx():
                    # Only include if the connecting bond is single.
                    next_bond = mol.GetBondBetweenAtoms(current_atom.GetIdx(), nbr.GetIdx())
                    if next_bond and next_bond.GetBondType() == Chem.BondType.SINGLE:
                        next_carbons.append(nbr)
            
            if len(next_carbons) > 1:
                # Branching detected; this is not an unbranched acyl chain.
                valid_chain = False
                break
            elif len(next_carbons) == 0:
                # End of the chain reached.
                break
            else:
                # Proceed down the chain
                prev_atom = current_atom
                current_atom = next_carbons[0]
        
        # Finally, check that the extracted chain count is exactly 14.
        if valid_chain and chain_count == 14:
            return True, ("Found an ester group with a saturated, unbranched 14-carbon acyl chain (tetradecanoate) " +
                           "derived from tetradecanoic acid.")
    
    return False, "No ester group with an unbranched, saturated 14-carbon acyl chain (tetradecanoate) was detected."

# Example usage:
# test_smiles = "CCCCCCCCCCCCCC(=O)OC"  # Represents methyl tetradecanoate
# result, reason = is_tetradecanoate_ester(test_smiles)
# print(result, reason)