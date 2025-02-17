"""
Classifies: CHEBI:25627 octadecadienoic acid
"""
"""
Classifies: octadecadienoic acid – a straight-chain C18 polyunsaturated fatty acid 
having exactly 2 non-aromatic C=C double bonds and a terminal carboxylic acid group.
"""

from rdkit import Chem

def is_octadecadienoic_acid(smiles: str):
    """
    Determines whether a molecule (given by its SMILES string)
    is an octadecadienoic acid. The definition enforced here is that
      • the molecule has exactly 18 carbon atoms in total,
      • the molecule contains a terminal carboxylic acid group (SMARTS: "[CX3](=O)[OX2H1]"),
      • the carbon atoms (ignoring heteroatoms) form a single connected, linear (unbranched) chain:
          meaning that in the induced graph on carbons from the molecule,
          the connected component (which must include the acid carbon) has exactly 18 atoms,
          with exactly two endpoints (degree 1) and every other carbon has degree 2,
      • exactly 2 non-aromatic C=C bonds appear along the main carbon chain.
      
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        (bool, str): Tuple.
                     True with a positive reason if it is octadecadienoic acid,
                     otherwise False with reason.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for carboxylic acid functionality using a SMARTS query.
    acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "Missing carboxylic acid functionality"
    
    # We assume the acid carbon (the one in C(=O)O) is the first atom of the first match.
    acid_carbon_idx = acid_matches[0][0]
    
    # Count all carbon atoms in the molecule.
    carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbon_atoms) != 18:
        return False, f"Expected 18 carbon atoms but found {len(carbon_atoms)}"
    
    # Build an undirected graph using only carbon atoms.
    # The nodes are atom indices for carbons. An edge exists if two carbons are bonded.
    carbon_graph = {}
    for atom in carbon_atoms:
        idx = atom.GetIdx()
        # Only add neighbors that are carbons.
        nbrs = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        carbon_graph[idx] = nbrs

    # Get the connected component (set of indices) containing the acid carbon.
    visited = set()
    stack = [acid_carbon_idx]
    while stack:
        current = stack.pop()
        if current in visited:
            continue
        visited.add(current)
        for nbr in carbon_graph.get(current, []):
            if nbr not in visited:
                stack.append(nbr)
    if len(visited) != 18:
        return False, ("Carbon network is not a single, straight 18-carbon chain. "
                       f"Found a connected carbon subgraph of size {len(visited)}")
    
    # Check for linearity: in a linear chain, exactly two carbons should have degree 1 (endpoints)
    # and every other carbon must have degree 2.
    endpoints = [idx for idx, nbrs in carbon_graph.items() if len(nbrs) == 1]
    if len(endpoints) != 2:
        return False, ("Carbon chain is branched or cyclic. "
                       f"Expected 2 terminal carbons (degree=1) but found {len(endpoints)}")
    for idx, nbrs in carbon_graph.items():
        if len(nbrs) > 2:
            return False, f"Carbon atom with index {idx} shows branching (degree >2 in carbon network)"
    
    # Now, find the unique simple path between the two endpoints.
    # We perform a simple DFS search from one endpoint to the other in the induced carbon graph.
    start, end = endpoints[0], endpoints[1]
    path = []
    found = []

    def dfs(current, target, current_path, visited_dfs):
        if current == target:
            found.extend(current_path)
            return True
        for nbr in carbon_graph[current]:
            if nbr in visited_dfs:
                continue
            if dfs(nbr, target, current_path + [nbr], visited_dfs | {nbr}):
                return True
        return False

    if not dfs(start, end, [start], {start}):
        return False, "Failed to find a linear path in the carbon network"
    path = found
    if len(path) != 18:
        return False, f"Main chain length is {len(path)} instead of 18 carbons"

    # Count non-aromatic double bonds along the path.
    double_bond_count = 0
    for i in range(len(path) - 1):
        bond = mol.GetBondBetweenAtoms(path[i], path[i+1])
        if bond is None:
            return False, f"Missing bond between carbons {path[i]} and {path[i+1]}"
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and (not bond.GetIsAromatic()):
            double_bond_count += 1

    if double_bond_count != 2:
        return False, (f"Found {double_bond_count} non-aromatic C=C bonds along the main chain; "
                       "exactly 2 are required")
    
    return True, ("Molecule is a straight-chain C18 fatty acid with exactly 2 non-aromatic C=C bonds "
                  "and a terminal carboxylic acid group")

# For testing (will run when the module is executed directly)
if __name__ == "__main__":
    # Sample test cases (the list includes examples from successful cases and errors).
    test_smiles = [
        # True positives:
        "CCCCCC\\C=C/C=C/[C@H](O)CCCCCCCC(O)=O",  # 9(R)-HODE
        "OC(=O)CCCCCCC\\C=C/C=C\\CCCCCC",           # 9Z,11Z-octadecadienoic acid
        "CC/C=C\\C/C=C\\CC(C(CCCCCCCC(O)=O)O)O",     # 9,10-DiHODE
        # False positives (should be rejected because of branching/substituents):
        "C(=C\\C/C=C\\CCCCCO)CCCCCCCC(=O)O",         # 18-hydroxylinoleic acid
        "OC(=O)CCCCC/C=C/C=C\\CCCCCCCC",             # 7-trans,9-cis-octadecadienoic acid
        # False negative example: wrong number of carbons
        "O(C(CCCCCCC(O)=O)/C=C/C(=O)CCCCCCCCC(O)=O",  # (11E)-13-hydroxy-10-oxo-11-octadecenoic acid (extra carbon)
    ]
    for sm in test_smiles:
        result, reason = is_octadecadienoic_acid(sm)
        print(f"SMILES: {sm}\nResult: {result}, Reason: {reason}\n")