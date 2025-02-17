"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: Fatty acid methyl ester

A fatty acid methyl ester is defined as the carboxylic ester obtained by the formal 
condensation of a fatty acid with methanol. Here we require that the molecule 
(i) contains a methyl ester group as defined by the SMARTS "[CX3](=O)O[CH3]",
(ii) the acyl portion (the group attached to the carbonyl carbon, excluding the ester oxygen) 
    forms a contiguous chain of carbon atoms of at least 3 carbons,
(iii) the chain must not contain aromatic atoms and must be “simple” (i.e. almost all substituents 
    are trivial methyl branches – if any non-carbon substituent is detected the molecule is rejected),
(iv) and the overall molecule is not decorated beyond expectation.
    
If these conditions are met, the function returns True with an explanation.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    
    The algorithm:
      1. Searches for the methyl ester group via the SMARTS "[CX3](=O)O[CH3]".
      2. For each found ester substructure, it identifies the carbonyl carbon and then 
         the acyl side (the neighbor that is not the carbonyl oxygen).
      3. It then uses a depth-first search (DFS) to find the longest contiguous path 
         of carbon atoms (only aliphatic, non‐aromatic) starting from that acyl atom.
      4. If the chain length is at least 3 carbons then the chain is scrutinized:
           • any substituent (branching off a chain atom, other than the one linking back to the carbonyl)
             is examined. Non‐carbon branches or branches longer than one atom are a sign of extra decoration.
      5. The overall heavy atom count of the molecule is compared to an expected minimum 
         (chain_length + 4) – if there is too much extra functionality the chain is rejected.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a fatty acid methyl ester, False otherwise.
        str: Explanation for the classification decision.
    """
    
    # Try to parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define the methyl ester SMARTS
    ester_smarts = "[CX3](=O)O[CH3]"
    ester_pattern = Chem.MolFromSmarts(ester_smarts)
    if ester_pattern is None:
        return False, "Error defining ester SMARTS"
    
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if not ester_matches:
        return False, "No methyl ester group found in the molecule"
    
    # DFS helper: from a given atom, follow bonds to atoms that are carbons, not in a ring, and non aromatic.
    def dfs_longest_chain(atom, visited):
        # Current path includes this atom.
        best_path = [atom.GetIdx()]
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() != 6:
                continue
            # We want only non-aromatic carbons (i.e. aliphatic chain) 
            # so if this neighbor is marked aromatic, skip it.
            if nbr.GetIsAromatic():
                continue
            if nbr.GetIdx() in visited:
                continue
            new_visited = visited | {nbr.GetIdx()}
            candidate = dfs_longest_chain(nbr, new_visited)
            candidate = [atom.GetIdx()] + candidate
            if len(candidate) > len(best_path):
                best_path = candidate
        return best_path
    
    # Examine each ester match.
    # The SMARTS match returns a tuple: (carbonylC, carbonylO, esterO, methylC)
    for match in ester_matches:
        carbonyl_idx, carbonylO_idx, esterO_idx, methyl_idx = match
        carbonyl_atom = mol.GetAtomWithIdx(carbonyl_idx)
        
        # Identify the acyl neighbor (attached to the carbonyl but not the carbonyl oxygen)
        acyl_neighbor = None
        for nbr in carbonyl_atom.GetNeighbors():
            if nbr.GetIdx() == carbonylO_idx:
                continue
            if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                acyl_neighbor = nbr
                break
        if acyl_neighbor is None:
            continue  # try next ester match
        
        # Compute the longest contiguous acyl chain starting from acyl_neighbor.
        visited = {carbonyl_idx}  # exclude the carbonyl carbon
        chain_path = dfs_longest_chain(acyl_neighbor, visited)
        chain_length = len(chain_path)
        
        if chain_length < 3:
            # A minimal fatty acyl fragment should have at least 3 carbons.
            continue
        
        # Check for extra decoration on the acyl chain.
        # For each atom in the chain, check neighbors that are not in the chain and not the carbonyl (for the first atom).
        extra_branch = False
        for i, atom_idx in enumerate(chain_path):
            atom = mol.GetAtomWithIdx(atom_idx)
            for nbr in atom.GetNeighbors():
                # Skip atoms that are part of the main chain.
                if nbr.GetIdx() in chain_path:
                    continue
                # Also allow the carbonyl atom if we are at the first atom.
                if i == 0 and nbr.GetIdx() == carbonyl_idx:
                    continue
                # We allow a trivial (methyl) branch if it is a carbon that itself has no further heavy neighbors.
                if nbr.GetAtomicNum() == 6 and not nbr.GetIsAromatic():
                    # Count heavy neighbors (excluding hydrogens)
                    heavy_nbrs = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1 and a.GetIdx() not in chain_path]
                    if len(heavy_nbrs) > 0:
                        extra_branch = True
                        break
                else:
                    # Any branch that is not carbon is not allowed.
                    extra_branch = True
                    break
            if extra_branch:
                break
        if extra_branch:
            # The chain shows decoration beyond a simple fatty acid
            continue
        
        # Check overall molecule simplicity.
        # Expected heavy atoms: acyl chain (chain_length) + carbonyl carbon (1) + carbonyl oxygen (1) 
        #                    + ester oxygen (1) + methyl carbon (1) = chain_length + 4.
        expected_heavy_atoms = chain_length + 4
        actual_heavy_atoms = mol.GetNumHeavyAtoms()
        # We allow a tolerance of +6 heavy atoms only.
        if actual_heavy_atoms > expected_heavy_atoms + 6:
            continue
        
        # Optionally compute the number of rotatable bonds in the acyl chain.
        rotatable_in_chain = 0
        for j in range(len(chain_path)-1):
            bond = mol.GetBondBetweenAtoms(chain_path[j], chain_path[j+1])
            if bond is None:
                continue
            if bond.GetBondTypeAsDouble() == 1.0 and not bond.IsInRing():
                rotatable_in_chain += 1
        
        overall_rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
        mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
        
        reason = (f"Found a methyl ester group with an acyl chain of {chain_length} carbon(s); "
                  f"rotatable bonds in chain: {rotatable_in_chain}, overall rotatable bonds: {overall_rot_bonds}, "
                  f"molecule heavy atoms: {actual_heavy_atoms} (expected ~{expected_heavy_atoms}), "
                  f"molecular weight: {mol_wt:.1f} Da")
        return True, reason
    
    return False, "Molecule contains a methyl ester group but the acyl portion does not meet fatty acid criteria"


# Example usage (testing several provided SMILES):
if __name__ == "__main__":
    test_smiles = [
        "CCCCCC1C(C/C(=C(\\C\\C(=C(\\[H])/CCCCCCC(=O)OC)\\[H])/[H])/[H])O1",  # 13-(3-pentyl-2-oxiranyl)-... (TP)
        "O[C@H](CCCCCCCC(OC)=O)/C=C/C=C/CCCCC",  # 9R-hydroxy-10E,12E-octadecadienoic acid methyl ester (TP)
        "O=C(OC)CCCCCCCCCCCC",  # methyl tridecanoate (TP)
        "O(C(=O)CCC(CCCC)C)C",  # formyl 4-methyl-octanoate (TP)
        "COC(=O)[C@@H]1O[C@@]1(C)CC\\C=C(/C)CC[C@H]1OC1(C)C",  # juvenile hormone III skipped bisepoxide (TP)
        "COC(=O)CCC(=O)CN",  # methyl 5-aminolevulinate (FP expected to be rejected)
        "COC(=O)\\C=C(/C)C(=O)OC",  # Dimethyl citraconate (FP expected to be rejected)
        "O=C(OC)CCC(=O)NCCCCCN(O)C(=O)CC(C)C"  # Tenacibactin A (FP expected to be rejected)
    ]
    
    for sm in test_smiles:
        valid, msg = is_fatty_acid_methyl_ester(sm)
        print(f"SMILES: {sm}\nResult: {valid}\nReason: {msg}\n")