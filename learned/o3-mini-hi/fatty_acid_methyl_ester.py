"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: Fatty acid methyl ester
A fatty acid methyl ester is defined as the carboxylic ester obtained by the formal 
condensation of a fatty acid (with its acyl chain) with methanol.
This algorithm first finds a methyl ester group via the SMARTS [CX3](=O)O[CH3]. 
Then for the ester match it identifies the acyl chain attached to the carbonyl carbon 
(excluding the ester oxygen) by finding the longest contiguous carbon-only chain 
(starting from that neighbor). It then computes:
  - chain_length (number of carbons in the acyl fragment)
  - internal rotatable bonds count (only bonds not at the chain ends)
Finally it checks that the overall molecule is “simple” (i.e. its heavy atom count is not 
excessively above the minimal expected count, which would indicate additional functionality).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    
    The algorithm:
      1. Searches for the methyl ester group using SMARTS "[CX3](=O)O[CH3]".
      2. For each match, identifies the acyl chain attached (the neighbor of the carbonyl 
         that is not the carbonyl oxygen). Then, it uses a depth-first search to select 
         the longest continuous carbon-only chain (allowing unsaturation).
      3. For chains with at least 5 carbons, it requires that the internal rotatable bonds 
         (bonds connecting non-terminal carbons) are at least half the chain length.
      4. It then checks that the overall molecule is not “decorated” with much extra functionality.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a fatty acid methyl ester, False otherwise.
        str: Explanation for the classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the SMARTS for a methyl ester group.
    ester_smarts = "[CX3](=O)O[CH3]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error defining ester SMARTS"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No methyl ester group found in the molecule"
    
    # Helper: DFS to find the longest simple path (list of atom indices) starting from a given atom.
    # We restrict the search to carbon atoms only.
    def dfs_longest_path(atom, visited):
        # current path includes this atom; start with length 1.
        best_path = [atom.GetIdx()]
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            if nbr.GetIdx() in visited:
                continue
            # Explore deeper: add current atom to visited.
            new_visited = visited.union({nbr.GetIdx()})
            candidate_path = dfs_longest_path(nbr, new_visited)
            candidate_path = [atom.GetIdx()] + candidate_path
            if len(candidate_path) > len(best_path):
                best_path = candidate_path
        return best_path

    # For each methyl ester match.
    # The SMARTS match returns (carbonylC, carbonylO, esterO, methylC)
    for match in ester_matches:
        carbonyl_idx, carbonylO_idx, esterO_idx, methyl_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the acyl neighbor of the carbonyl (exclude the carbonyl oxygen).
        acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == carbonylO_idx:
                continue
            if nbr.GetAtomicNum() == 6:
                acyl_neighbor = nbr
                break
        if acyl_neighbor is None:
            continue  # Try next ester match if no acyl chain found.
        
        # Compute the longest linear chain (as a list of atom indices) starting from the acyl neighbor.
        # We ensure not to go back to the carbonyl.
        visited = {carbonyl_idx}  # exclude the carbonyl from the chain.
        longest_chain = dfs_longest_path(acyl_neighbor, visited)
        chain_length = len(longest_chain)
        
        if chain_length < 3:
            # Even a minimal fatty acid should have at least three carbons in its acyl chain.
            continue
        
        # Count internal rotatable bonds along the longest chain.
        # For a linear chain, a bond is considered internal if it connects two non‐terminal atoms.
        # We count a bond if it is a single bond and not in a ring.
        rotatable_in_chain = 0
        if chain_length >= 3:
            # In a linear chain, internal bonds are those not incident to the endpoints.
            for i in range(1, chain_length-1):
                # Get bond connecting the atoms at positions i and i+1 in the chain path.
                a_idx = longest_chain[i]
                b_idx = longest_chain[i+1]
                bond = mol.GetBondBetweenAtoms(a_idx, b_idx)
                if bond is None:
                    continue
                if bond.GetBondTypeAsDouble() != 1.0:
                    continue
                if bond.IsInRing():
                    continue
                # In a simple chain the internal atoms (i not 0 and not chain_length-1) are not terminal.
                rotatable_in_chain += 1

        # For chains that are long (>=5 carbons), we require a minimum number of rotatable bonds.
        if chain_length >= 5:
            required_rot_bonds = chain_length // 2  # heuristic: at least half of carbons (excluding terminals)
            if rotatable_in_chain < required_rot_bonds:
                # Possibly the chain is too rigid (e.g. if many bonds are double bonds or constrained)
                continue
        
        # Next, check for the "simplicity" of the molecule.
        # For a simple fatty acid methyl ester, the expected heavy atom count should be roughly:
        #   acyl chain (chain_length) + carbonyl carbon (1) + ester O (1) + carbonyl oxygen (1) + methyl carbon (1) = chain_length + 4.
        expected_heavy_atoms = chain_length + 4
        actual_heavy_atoms = mol.GetNumHeavyAtoms()
        # If the molecule has much more than expected (say >150%), then there is extra decoration.
        if actual_heavy_atoms > expected_heavy_atoms * 1.5:
            continue
        
        # Also check the overall number of rotatable bonds for reasonability.
        overall_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        
        reason = (f"Found a methyl ester group with an acyl chain of {chain_length} carbon(s); "
                  f"internal rotatable bonds in chain: {rotatable_in_chain}, overall rotatable bonds: {overall_rot_bonds}, "
                  f"molecule heavy atoms: {actual_heavy_atoms} (expected ~{expected_heavy_atoms}), "
                  f"molecular weight: {mol_wt:.1f} Da")
        return True, reason
    
    return False, "Molecule contains a methyl ester group but the acyl portion does not meet fatty acid criteria"

# Example usage:
if __name__ == "__main__":
    # Some test SMILES (the list below includes a few of the examples provided)
    test_smiles = [
        "CCCCCC1C(C/C(=C(\\C\\C(=C(\\[H])/CCCCCCC(=O)OC)\\[H])/[H])/[H])O1",  # 13-(3-pentyl-2-oxiranyl)-...
        "O[C@H](CCCCCCCC(OC)=O)/C=C/C=C/CCCCC",  # 9R-hydroxy-10E,12E-octadecadienoic acid, methyl ester
        "O=C(OC)CCCCCCCCCCCC",  # methyl tridecanoate
        "COC(=O)[C@@H]1O[C@@]1(C)CC\\C=C(/C)CC[C@H]1OC1(C)C",  # juvenile hormone III skipped bisepoxide (was false negative before)
        "COC(=O)CCC(=O)CN"  # methyl 5-aminolevulinate (false positive in previous attempt)
    ]
    for sm in test_smiles:
        valid, msg = is_fatty_acid_methyl_ester(sm)
        print(f"SMILES: {sm}\nResult: {valid}\nReason: {msg}\n")