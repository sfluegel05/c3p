"""
Classifies: CHEBI:27208 unsaturated fatty acid
"""
"""
Classifies: Unsaturated Fatty Acid
Definition (heuristic):
  A molecule is deemed an unsaturated fatty acid if it meets all of the following:
    1. The SMILES string represents a valid molecule.
    2. It contains exactly one carboxylic acid group – as recognized by the SMARTS "C(=O)[O;H1,O-]".
    3. The carboxyl (acid) carbon is terminal; that is, it is attached to exactly one carbon atom.
    4. Starting from that attached carbon, if one finds the longest strictly carbon-based simple path,
       its length (number of carbons) is at least MIN_CHAIN (here we use 4 as a relaxed requirement).
    5. Along that chain at least one bond is a double bond (C=C).
    6. The molecule does not contain any carbon–carbon triple bonds.
    
Note:
  – This is a heuristic approach. Some fatty acids (e.g. with rings or unusual branching) may be mis‐classified.
  – The code “walks” from the acid’s neighbor along bonds between carbon atoms only.
  
Examples:
  • resolvin D6, leukotriene A4, farnesoic acid ... are accepted
  • Molecules with extra acid groups, very short chains (e.g. isocrotonic acid), or triple bonds are rejected.
  
Improvements over the previous approach:
  – The longest chain is computed explicitly “walking” and the bond orders along that path are examined.
  – The chain length threshold is relaxed so that known short-chain unsaturated acids (e.g. pent‐4‐enoic acid)
    are successfully classified, while still eliminating many false positives.
"""

from rdkit import Chem

def is_unsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is an unsaturated fatty acid based on its SMILES string.
    
    The criteria are:
      1. Molecule must be valid.
      2. Must contain exactly one carboxylic acid group defined by SMARTS "C(=O)[O;H1,O-]".
      3. The carboxyl carbon (from the acid group) must be terminal – attached to exactly one carbon.
      4. The longest simple carbon-only path (starting from the acid’s only carbon neighbor)
         must be at least MIN_CHAIN carbons long.
      5. At least one bond in that path must be a C=C double bond.
      6. The molecule must not contain any carbon–carbon triple bonds.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule qualifies as an unsaturated fatty acid, False otherwise.
        str: Explanation of the result.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Reject molecules with any carbon–carbon triple bonds.
    triple_bond_pattern = Chem.MolFromSmarts("C#C")
    if mol.HasSubstructMatch(triple_bond_pattern):
        return False, "Contains carbon–carbon triple bonds which are not allowed"
    
    # Define SMARTS for a carboxylic acid group (either protonated or deprotonated)
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H1,O-]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Does not contain a carboxylic acid group"
    if len(acid_matches) != 1:
        return False, "Molecule should have exactly one carboxylic acid group"
    
    # Our SMARTS has the acid carbon as atom 0.
    acid_match = acid_matches[0]
    acid_c_idx = acid_match[0]
    acid_c = mol.GetAtomWithIdx(acid_c_idx)
    
    # For a terminal acid group, the acid carbon must be attached to exactly one carbon.
    carbon_neighbors = [nbr.GetIdx() for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxyl group is not terminal (acid carbon is attached to more than one carbon)"
    
    chain_start_idx = carbon_neighbors[0]
    
    # Define a DFS routine to find the longest carbon-only simple path starting from a given atom.
    # We also collect the bond orders along the path.
    def dfs(current_idx, visited):
        visited.add(current_idx)
        current_atom = mol.GetAtomWithIdx(current_idx)
        # Store as (max_length, bond_info, path) where bond_info is a list of bond orders traversed.
        best = (1, [], [current_idx])
        for nbr in current_atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            nbr_idx = nbr.GetIdx()
            if nbr_idx in visited:
                continue
            # Get the bond between current and neighbor.
            bond = mol.GetBondBetweenAtoms(current_idx, nbr_idx)
            bond_order = bond.GetBondTypeAsDouble()  # 1.0 for single, 2.0 for double, etc.
            new_length, new_bonds, new_path = dfs(nbr_idx, visited.copy())
            candidate_length = 1 + new_length
            candidate_bonds = [bond_order] + new_bonds
            candidate_path = [current_idx] + new_path
            if candidate_length > best[0]:
                best = (candidate_length, candidate_bonds, candidate_path)
        return best

    # Set a minimum chain length (number of carbons in the chain, not counting the acid carbon).
    MIN_CHAIN = 4
    chain_length, bonds_along_chain, path = dfs(chain_start_idx, set())
    
    if chain_length < MIN_CHAIN:
        return False, f"Alkyl chain too short for a fatty acid (chain length = {chain_length})"
    
    # Check for unsaturation on the main chain: At least one bond in the chain should be a double (order == 2)
    if not any(bond_order == 2.0 for bond_order in bonds_along_chain):
        return False, "Main carbon chain does not contain any C=C double bonds (unsaturation missing)"
    
    # All criteria met – accept the molecule.
    return True, (f"Contains a terminal carboxyl group, a sufficiently long carbon chain (length = {chain_length}), "
                  "and at least one C=C double bond in the main chain (unsaturation).")

# Example usage (for testing purposes):
if __name__ == '__main__':
    test_smiles = [
        # True positives:
        "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O",  # resolvin D6
        "CCCC\\C=C/C\\C=C/C=C/C=C/[C@@H]1O[C@H]1CCCC(O)=O",  # leukotriene A4
        "CCCCCCC/C=C/CCCCCCC/C(=O)O",  # (2E,4E,6E,8E,10E)-octadecapentaenoic acid-like
        "CC\\C(CC[C@H]1O[C@@]1(C)CC)=C/CC\\C(C)=C\\C(O)=O",  # juvenile hormone I acid
        "OC(=O)CCC=C",  # pent-4-enoic acid (chain length 4; now accepted)
        "[H]C(C)=CC([H])=CC(O)=O",  # sorbic acid (chain length 3 or 4 depending on interpretation)
        # False positive example:
        "CCCCCCC\\C=C\\CCCCCCCCC([O-])=O",  # (10E)-octadecenoate – may be a salt derivative and complex
    ]
    
    for smi in test_smiles:
        result, reason = is_unsaturated_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n{'-'*50}")