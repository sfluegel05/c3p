"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: Octanoate Ester 
Definition: Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).
In other words, the ester contains an acyl group that is CH3(CH2)6C(=O)–.
This program looks for each ester (i.e. –C(=O)O–) in the molecule and checks whether 
its acyl part is an unbranched chain of 7 carbons (CH2 repeated 6 times and ending in a CH3),
so that when counting the carbonyl carbon the chain is 8 carbons total.
"""

from rdkit import Chem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.
    
    An octanoate ester is defined as any fatty acid ester in which the acyl part is 
    CH3(CH2)6C(=O)–. In an ideal octanoate ester the ester bond is formed via a carboxylic
    acid (octanoic acid) whose acyl chain (excluding the carbonyl) is exactly 7 carbons long:
    the first 6 should be methylene groups (–CH2–) and the terminal one a methyl group (–CH3).
    
    To improve upon the previous strategy we first find all ester groups (using the SMARTS "C(=O)O")
    and then for each ester we check that:
      1. The carbonyl carbon (C(=O)) has exactly one carbon neighbor (acyl side) besides the carbonyl oxygen.
      2. Following that acyl chain gives exactly 7 carbons,
         with the first 6 carbons being CH2 and the terminal carbon a CH3.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if all ester groups (if any) are derived from octanoic acid, False otherwise.
        str: A reason for the classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all ester groups using a SMARTS pattern.
    # This pattern matches a carbonyl (C(=O)) directly bonded to an oxygen.
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern, uniquify=True)
    
    # If no ester bond is present then by definition the molecule is not an ester.
    if not ester_matches:
        return False, "No ester group found"
    
    # Helper: check if a given carbon atom is CH2 (2 hydrogens) or CH3 (3 hydrogens)
    def is_ch2(atom):
        # Only count explicit + implicit hydrogens.
        return atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 2

    def is_ch3(atom):
        return atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() == 3
    
    # For each ester group found, we will verify the acyl chain.
    for match in ester_matches:
        # match[0] = carbonyl C; match[1] = carbonyl O (double-bonded); match[2] = ester O (single-bonded)
        carbonyl_atom = mol.GetAtomWithIdx(match[0])
        ester_oxygen = mol.GetAtomWithIdx(match[2])
        
        # Find the acyl side: among neighbors of carbonyl, skip the carbonyl oxygen (match[1]) and also skip the ester oxygen.
        acyl_neighbors = []
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() not in (match[1], match[2]) and nbr.GetAtomicNum() == 6:
                acyl_neighbors.append(nbr)
        if len(acyl_neighbors) != 1:
            return False, "Ester group found with ambiguous acyl connectivity"
        acyl_atom = acyl_neighbors[0]
        
        # Now traverse the acyl chain from acyl_atom.
        # For an octanoate chain derived from octanoic acid,
        # the chain (excluding the carbonyl atom) should have exactly 7 carbons:
        # positions 1-6: CH2 groups and position7: CH3 (terminal) [all in a linear, unbranched chain]
        chain_atoms = []
        current_atom = acyl_atom
        previous_atom = carbonyl_atom
        while True:
            chain_atoms.append(current_atom)
            # Look for the next carbon atom that is connected linearly (exclude the atom we came from).
            next_carbons = [nbr for nbr in current_atom.GetNeighbors() 
                            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() != previous_atom.GetIdx()]
            # For an unbranched chain, expect at most one next carbon.
            if len(next_carbons) == 0:
                # End of chain reached.
                break
            elif len(next_carbons) > 1:
                return False, "Branched acyl chain found"
            else:
                previous_atom, current_atom = current_atom, next_carbons[0]
        
        # We require exactly 7 carbons in the acyl chain.
        if len(chain_atoms) != 7:
            return False, f"Ester group acyl chain length is {len(chain_atoms)} (expected 7)"
        
        # Check that the first 6 atoms are CH2 and the last is CH3.
        for i, atom in enumerate(chain_atoms):
            if i < 6:
                if not is_ch2(atom):
                    return False, "Acyl chain does not have six consecutive CH2 groups"
            else:  # last atom
                if not is_ch3(atom):
                    return False, "Acyl chain terminal group is not CH3"
    
    # If we got here, every ester group found is an octanoate ester.
    return True, "All ester groups are derived from octanoic acid (octanoate ester)."

# Example usage:
if __name__ == '__main__':
    test_smiles = [
        "CCCCCCCC(=O)OC[C@H](O)CO",  # 3-octanoyl-sn-glycerol (true positive)
        "CCCCCCCC(=O)OC",            # methyl octanoate (true positive)
        "O(C(=O)CCCCCCC)C(C)C",       # isopropyl octanoate (may be flagged as false because of branching in the alcohol part is allowed, but our check focuses only on the acyl chain)
    ]
    for sm in test_smiles:
        flag, reason = is_octanoate_ester(sm)
        print(f"SMILES: {sm}\nClassification: {flag}\nReason: {reason}\n")