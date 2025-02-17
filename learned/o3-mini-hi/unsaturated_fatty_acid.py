"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
#!/usr/bin/env python
"""
Classifies: Unsaturated Fatty Acid
Definition (heuristic):
  A molecule is considered an unsaturated fatty acid if:
    1. It is a valid molecule.
    2. It contains exactly one free (protonated) carboxylic acid group.
         (Here we match a carboxyl group using the SMARTS "C(=O)[OH]".)
    3. The acid carbon (the carbon of the acid group) is terminal (i.e. it is connected to exactly one carbon).
    4. Starting from that single carbon neighbor, the longest simple path of carbons is at least a minimum length.
    5. At least one bond along that carbon chain is a C=C double bond.
    6. The molecule does not contain any carbon–carbon triple bonds.
    
This approach was refined to address previous false positives (e.g. deprotonated salts were accepted)
and false negatives (molecules with extra acid groups were flagged) by enforcing both:
  – That the acid group is protonated (free acid) and unique.
  – That the chain length and unsaturation are computed from the carbon attached to the acid carbon.
  
Note: This is a heuristic method; some fatty acids with highly branched chains or additional groups 
may not be handled perfectly.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    
    The criteria are:
      1. The molecule can be parsed.
      2. It must not contain any carbon–carbon triple bonds.
      3. It must contain exactly one free (protonated) carboxylic acid group
         as defined by the SMARTS "C(=O)[OH]".
      4. The carboxyl carbon must be terminal (i.e. attached to exactly one other carbon).
      5. The longest simple carbon-only path starting from that neighbor must have at least MIN_CHAIN carbons.
      6. At least one bond along that chain must be a double (C=C) bond.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as an unsaturated fatty acid, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with any carbon–carbon triple bonds.
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon–carbon triple bonds which are not allowed"
    
    # Define a SMARTS for a free (protonated) carboxylic acid group.
    # This pattern will match a carbonyl with an -OH (not O-).
    acid_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Does not contain a free (protonated) carboxylic acid group"
    if len(acid_matches) != 1:
        return False, "Molecule should have exactly one free (protonated) carboxylic acid group"
    
    # In the SMARTS, we assume atom 0 is the acid carbon.
    acid_match = acid_matches[0]
    acid_c_idx = acid_match[0]
    acid_c = mol.GetAtomWithIdx(acid_c_idx)
    
    # The acid carbon must be terminal.  Count its carbon neighbors.
    carbon_neighbors = [nbr.GetIdx() for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxyl group is not terminal (acid carbon is attached to more than one carbon)"
    
    # Starting from the single carbon neighbor, compute the longest simple carbon-only chain.
    chain_start_idx = carbon_neighbors[0]
    
    # DFS routine: returns a tuple (chain_length, bond_orders, path)
    # chain_length counts the number of carbons visited (starting at the current atom = 1)
    def dfs(current_idx, visited):
        visited.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        best = (1, [], [current_idx])  # (length, list of bond orders taken, atom path)
        # Explore neighboring carbons not already visited.
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            bond = mol.GetBondBetweenAtoms(current_idx, nbr_idx)
            bond_order = bond.GetBondTypeAsDouble()  # 1.0 for single, 2.0 for double, etc.
            new_length, new_bond_orders, new_path = dfs(nbr_idx, visited.copy())
            candidate_length = 1 + new_length
            candidate_bond_orders = [bond_order] + new_bond_orders
            candidate_path = [current_idx] + new_path
            if candidate_length > best[0]:
                best = (candidate_length, candidate_bond_orders, candidate_path)
        return best
    
    # Set a minimum required chain length (number of carbons in the chain, excluding the acid carbon).
    MIN_CHAIN = 4
    chain_length, bonds_along_chain, path = dfs(chain_start_idx, set())
    
    if chain_length < MIN_CHAIN:
        return False, f"Alkyl chain too short for a fatty acid (chain length = {chain_length})"
    
    # Check for unsaturation: At least one bond in the chain must be a double bond.
    if not any(bond_order == 2.0 for bond_order in bonds_along_chain):
        return False, "Main carbon chain does not contain any C=C double bonds (unsaturation missing)"
    
    return True, (f"Contains a terminal carboxyl group, a sufficiently long carbon chain (length = {chain_length}), "
                  "and at least one C=C double bond in the main chain (unsaturation).")

# Example usage (for testing purposes only):
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O",  # resolvin D6
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CCCCCCC/C=C/CCCCCCC/C(=O)O",  # (2E,4E,6E,8E,10E)-octadecapentaenoic acid-like
        "CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O",  # juvenile hormone I acid
        "C(C(O)=O)C/C=C\\C[C@H](\\C=C\\C=C/C=C/C=C/[C@H]([C@@H](C/C=C\\CC)O)O)O",  # aspirin-triggered resolvin D2
        "CC(C)=CCCC(C)=CCCC(C)=CC(O)=O",  # farnesoic acid
        "OC(=O)CCC=C",  # pent-4-enoic acid (chain length 4; now accepted)
        "[H]C(C)=CC([H])=CC(O)=O",  # sorbic acid (if interpreted as chain length 4)
        # False positive example (should be rejected due to being a deprotonated salt):
        "CCCCCCC\\C=C\\CCCCCCCCC([O-])=O",  # (10E)-octadecenoate
    ]
    
    for smi in test_smiles:
        result, reason = is_unsaturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*50}")