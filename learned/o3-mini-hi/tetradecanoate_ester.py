"""
Classifies: CHEBI:87691 tetradecanoate ester
"""
#!/usr/bin/env python
"""
Classifies: Tetradecanoate ester
A tetradecanoate ester is defined as a fatty acid ester formed by condensation of the carboxyl group of tetradecanoic acid 
(with a 14‐carbon (myristoyl) acyl group, CH3(CH2)12CO–) with a hydroxy group of an alcohol or phenol.
"""

from rdkit import Chem

def is_tetradecanoate_ester(smiles: str):
    """
    Determines if a molecule is a tetradecanoate ester based on its SMILES string.
    
    A tetradecanoate ester has at least one ester group (R–C(=O)–O–R')
    where the acyl (fatty acid) portion corresponds to tetradecanoic acid 
    (i.e. an unbranched chain of 14 carbons, with the carbonyl carbon counting as one).
    
    The function searches for a carbonyl carbon that is part of an ester group:
       - The carbonyl carbon should be bound via a double bond to O (carbonyl oxygen)
         and via a single bond to an oxygen (ester oxygen).
       - The remaining neighbor (should be a carbon) is assumed to be the start of the acyl chain.
       - We then traverse the acyl chain linearly, ensuring that each carbon (except the terminus)
         connects to exactly one additional carbon. If the total number of carbons (including the carbonyl) is 14,
         we return True.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        (bool, str): Tuple with True and a reason if a tetradecanoate ester group is found, 
                     or False with a reason if not.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Iterate through all atoms to look for a candidate ester carbonyl carbon.
    for atom in mol.GetAtoms():
        # Only consider carbon atoms
        if atom.GetAtomicNum() != 6:
            continue

        bonds = atom.GetBonds()
        carbonyl_oxygen = None
        ester_oxygen = None
        acyl_neighbor = None
        
        # We expect that in an ester group, the carbonyl carbon is attached to:
        #   - one oxygen by a DOUBLE bond (the carbonyl oxygen)
        #   - one oxygen by a SINGLE bond (the ester oxygen)
        #   - one carbon atom (the start of the acyl chain)
        for bond in bonds:
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 8:
                if bond.GetBondType() == Chem.BondType.DOUBLE:
                    carbonyl_oxygen = nbr
                elif bond.GetBondType() == Chem.BondType.SINGLE:
                    ester_oxygen = nbr
            elif nbr.GetAtomicNum() == 6:
                acyl_neighbor = nbr
        
        # Check if all three parts are present
        if not (carbonyl_oxygen and ester_oxygen and acyl_neighbor):
            continue
        
        # We have an ester carbonyl candidate.
        # Next, traverse the suspected acyl chain starting from the carbonyl carbon.
        # Our expected acyl chain should be linear: CH3-(CH2)12-CO–
        # Count the carbonyl carbon as part of the chain.
        chain_count = 1  # start counting the carbonyl carbon itself
        prev_atom = atom
        current_atom = acyl_neighbor
        valid_chain = True
        
        # Traverse linearly
        while True:
            # Check that current_atom is carbon
            if current_atom.GetAtomicNum() != 6:
                valid_chain = False
                break

            chain_count += 1
            
            # Determine the next carbon in the chain.
            # Get all carbon neighbors of current_atom (ignoring hydrogens) that are not the one we came from.
            next_carbons = [nbr for nbr in current_atom.GetNeighbors() 
                            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != prev_atom.GetIdx()]
            
            # For a linear unbranched chain, the intermediate carbon should have exactly one carbon neighbor (the next one).
            if len(next_carbons) > 1:
                # The chain is branched; not a simple fatty acyl chain.
                valid_chain = False
                break
            elif len(next_carbons) == 0:
                # We've reached a terminal carbon.
                break
            else:
                # Proceed to next carbon in the chain.
                prev_atom = current_atom
                current_atom = next_carbons[0]
        
        # For tetradecanoate the total number of carbons in the acyl chain (including the carbonyl carbon) must be 14.
        if valid_chain and chain_count == 14:
            return True, "Found an ester group with a 14-carbon acyl chain (tetradecanoate) derived from tetradecanoic acid."
    
    return False, "No ester group with a 14-carbon (tetradecanoate) acyl chain was detected."

# Example usage (uncomment for testing):
# test_smiles = "CCCCCCCCCCCCCC(=O)OC"  # This represents methyl tetradecanoate
# result, reason = is_tetradecanoate_ester(test_smiles)
# print(result, reason)