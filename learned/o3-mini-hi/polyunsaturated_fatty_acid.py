"""
Classifies: CHEBI:26208 polyunsaturated fatty acid
"""
"""
Classifies: Polyunsaturated Fatty Acid
Definition: A fatty acid with a terminal –COOH group that is attached to a single, 
            linear, non‐aromatic aliphatic chain and containing more than one 
            non‐aromatic carbon–carbon double bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_polyunsaturated_fatty_acid(smiles: str):
    """
    Determines if a molecule is a polyunsaturated fatty acid based on its SMILES string.
    The molecule must contain a terminal carboxylic acid group (the acid carbon is
    attached to exactly one carbon) and a long, linear, non-aromatic fatty acyl chain 
    with more than one non-aromatic C=C double bond.
    
    We additionally:
      - Extract the linear chain using a DFS that only follows non-aromatic, non-ring carbon atoms.
      - Only accept the chain if (a) it has a minimum length (here 8 carbons) and (b) it is strictly linear (each “internal” atom in the chain is connected only to two chain atoms).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule qualifies as a polyunsaturated fatty acid, False otherwise.
        str: Explanation for classification decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Look for a carboxylic acid group using SMARTS.
    # This pattern matches a carbon with two oxygens (one as =O, one as -OH).
    ca_smarts = "[CX3](=O)[O;H]"
    ca_group = Chem.MolFromSmarts(ca_smarts)
    ca_matches = mol.GetSubstructMatches(ca_group)
    if not ca_matches:
        return False, "No carboxylic acid group found; not a fatty acid"
    
    # 2. Check for a terminal acid: choose one where the acid carbon (the C in the COOH, index 0)
    # is attached to exactly one carbon (i.e. the chain-start).
    terminal_acid_found = False
    acid_atom = None
    chain_start = None
    for match in ca_matches:
        # In our match, the first atom is the carbonyl carbon.
        acid_c = mol.GetAtomWithIdx(match[0])
        carbon_neighbors = [nbr for nbr in acid_c.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            acid_atom = acid_c
            chain_start = carbon_neighbors[0]
            break
    if not terminal_acid_found or chain_start is None:
        return False, "Carboxylic acid group not terminal; not a typical fatty acid"
    
    # 3. Extract the chain. We want a linear chain of carbon atoms that are:
    #    - non-aromatic and not in any ring.
    # We perform a DFS that finds all simple paths starting from chain_start.
    # We then choose the longest chain that satisfies that all internal atoms are linked to exactly 2 chain atoms.
    def dfs_chain(atom, coming_from, path):
        """
        Recursive DFS: extend path with neighbors that are carbon (atomic num 6),
        non-aromatic, not in a ring, and not visited.
        """
        best = path
        for nbr in atom.GetNeighbors():
            if coming_from is not None and nbr.GetIdx() == coming_from.GetIdx():
                continue
            if nbr.GetAtomicNum() != 6 or nbr.GetIsAromatic() or nbr.IsInRing():
                continue
            if nbr.GetIdx() in path:
                continue
            candidate = dfs_chain(nbr, atom, path + [nbr.GetIdx()])
            if len(candidate) > len(best):
                best = candidate
        return best

    # Start from chain_start. We use a set of indices to avoid cycles.
    initial_path = [chain_start.GetIdx()]
    chain_path = dfs_chain(chain_start, acid_atom, initial_path)
    chain_length = len(chain_path)
    
    MIN_CHAIN_LENGTH = 8
    if chain_length < MIN_CHAIN_LENGTH:
        return False, f"Fatty acid chain length only {chain_length} carbons; too short to be typical"
    
    # 4. Check the linearity: in a linear chain, the two terminal carbons should have only one neighbor 
    # from the chain and all internal ones exactly two. We recreate the subgraph connectivity for the chain.
    chain_atoms = {idx: mol.GetAtomWithIdx(idx) for idx in chain_path}
    # Build neighbor count within the chain for each atom.
    neighbor_counts = {idx: 0 for idx in chain_path}
    for idx in chain_path:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetIdx() in chain_atoms:
                neighbor_counts[idx] += 1
    # Check: endpoints should have count==1; internal > endpoints should be ==2.
    endpoints = [idx for idx, count in neighbor_counts.items() if count == 1]
    internals = [idx for idx, count in neighbor_counts.items() if count != 1]
    if len(endpoints) != 2:
        # if not exactly two endpoints it means the chain is branched.
        return False, "Acyl chain is branched; expected linear chain for a typical fatty acid"
    for idx in internals:
        if neighbor_counts[idx] != 2:
            return False, "Acyl chain is not linear; found branching in the fatty chain"
    
    # 5. Count the number of non-aromatic double bonds along the chain.
    # We examine bonds between consecutive atoms in the chain as they appear along the longest path.
    # Note: The DFS path might not reflect the “natural” order but for a linear chain the bonds present in the path
    # should be the ones of interest.
    double_bond_count = 0
    for i in range(len(chain_path)-1):
        bond = mol.GetBondBetweenAtoms(chain_path[i], chain_path[i+1])
        if bond is not None:
            if bond.GetBondType() == Chem.BondType.DOUBLE and not bond.GetIsAromatic():
                double_bond_count += 1

    if double_bond_count <= 1:
        return False, f"Found {double_bond_count} non-aromatic C=C bond(s) in chain; need more than one to qualify as polyunsaturated"
    
    return True, f"Contains a terminal carboxylic acid group attached to a linear chain of {chain_length} carbons and {double_bond_count} non-aromatic double bonds; qualifies as a polyunsaturated fatty acid"


# Example usage:
if __name__ == "__main__":
    # Test one of the examples, resolvin D6.
    test_smiles = "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O"
    result, reason = is_polyunsaturated_fatty_acid(test_smiles)
    print(result, reason)