"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
#!/usr/bin/env python
"""
Classifies: Tetradecanoate ester

A tetradecanoate ester is defined as a fatty acid ester obtained by condensation of the carboxyl group 
of tetradecanoic acid (myristic acid; CH3(CH2)12CO–) with a hydroxy group of an alcohol or phenol.
Thus, the ester group must contain a saturated linear acyl chain of 14 carbons (with the carbonyl carbon included).
This program uses RDKit to search for an ester carbonyl that meets the following criteria:
  • The carbonyl carbon is bonded to exactly two oxygen atoms:
      – one via a double bond (the carbonyl oxygen)
      – one via a single bond (the ester oxygen that connects to the residue from the alcohol)
  • The carbonyl carbon is also bonded to a single carbon atom; this carbon is taken as the first atom
    of the acyl chain.
  • Traversing the acyl chain using only single bonds (and enforcing a linear, unbranched route),
    the total number of carbons (starting with the carbonyl carbon) must be exactly 14.
  • Any deviation (branching, unsaturation in the acyl chain, or length not equal to 14) causes this
    ester group to be rejected.
    
If any ester group in the molecule meets these criteria, the function returns True with a matching reason.
Otherwise, the molecule is not classified as a tetradecanoate ester.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    
    Args:
        smiles (str): SMILES representation of the molecule.
        
    Returns:
        (bool, str): Tuple where the first element is True if a tetradecanoate ester is found,
                     and the second element is a reason message. If not found, returns False and a message.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Loop over all atoms looking for candidate carbonyl carbons that could be part of an ester.
    for atom in mol.GetAtoms():
        # Only consider carbon atoms for the carbonyl carbon candidate.
        if atom.GetAtomicNum() != 6:
            continue
        
        bonds = atom.GetBonds()
        # Expect exactly two bonds from the candidate to oxygen atoms:
        #   • one double-bonded oxygen (carbonyl oxygen)
        #   • one single-bonded oxygen (ester oxygen)
        carbonyl_oxygen = None
        ester_oxygen = None
        acyl_neighbor = None  # the carbon that begins the fatty acyl chain
        
        # Loop over bonds attached to the candidate carbonyl carbon.
        for bond in bonds:
            nbr = bond.GetOtherAtom(atom)
            # If neighbor is oxygen, check bond type.
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    # Should be the carbonyl oxygen.
                    if carbonyl_oxygen is None:
                        carbonyl_oxygen = nbr
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    # Should be the ester oxygen.
                    if ester_oxygen is None:
                        ester_oxygen = nbr
            # A typical carbonyl carbon (in an ester) is also bonded to one carbon (the acyl chain start).
            elif nbr.GetAtomicNum() == 6:
                # We expect exactly one acyl-chain carbon
                if acyl_neighbor is None:
                    acyl_neighbor = nbr
                else:
                    # More than one carbon neighbor? Not the simple ester we want.
                    acyl_neighbor = None
                    break
        
        # Check that the candidate carbonyl meets the criteria.
        if not (carbonyl_oxygen and ester_oxygen and acyl_neighbor):
            continue  # move on to the next candidate
        
        # Now traverse the acyl chain.
        # The chain should be saturated, unbranched and should include exactly 14 carbons,
        # with the carbonyl carbon counted as the first.
        chain_count = 1  # start with the carbonyl carbon (counted as part of the acyl chain)
        prev_atom = atom
        current_atom = acyl_neighbor
        
        valid_chain = True
        # Traverse along single bonds only.
        while True:
            # Ensure the current atom is carbon and the connecting bond from prev_atom is single.
            bond = mol.GetBondBetweenAtoms(prev_atom.GetIdx(), current_atom.GetIdx())
            if bond is None or bond.GetBondType() != Chem.BondType.SINGLE:
                valid_chain = False
                break

            # Check that current_atom is sp3 (i.e. saturated) or at least not part of a double bond.
            # (We simply check that it does not have any double bonds to oxygen or other atoms in this chain.)
            # For our purposes the absence of multiple bonds along the acyl route is sufficient.
            if current_atom.GetAtomicNum() != 6:
                valid_chain = False
                break
            
            chain_count += 1
            
            # Get the neighbor carbons (exclude the one we just came from).
            next_carbons = []
            for nbr in current_atom.GetNeighbors():
                if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx():
                    # Only consider bonds that are SINGLE.
                    bnd = mol.GetBondBetweenAtoms(current_atom.GetIdx(), nbr.GetIdx())
                    if bnd and bnd.GetBondType() == Chem.BondType.SINGLE:
                        next_carbons.append(nbr)
            
            if len(next_carbons) > 1:
                # Branching detected; not a simple acyl chain.
                valid_chain = False
                break
            elif len(next_carbons) == 0:
                # End of the chain reached.
                break
            else:
                # Continue down the chain.
                prev_atom = current_atom
                current_atom = next_carbons[0]
        
        if valid_chain and chain_count == 14:
            return True, "Found an ester group with a saturated 14-carbon acyl chain (tetradecanoate) derived from tetradecanoic acid."
    
    return False, "No ester group with a saturated 14-carbon acyl chain (tetradecanoate) was detected."

# Example usage:
# test_smiles = "CCCCCCCCCCCCCC(=O)OC"  # Represents methyl tetradecanoate
# result, reason = is_tetradecanoate_ester(test_smiles)
# print(result, reason)