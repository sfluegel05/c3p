"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: Omega-hydroxy fatty acid
Definition: A naturally-occurring straight-chain fatty acid (acyclic, unbranched)
            composed only of C, H, and O that contains one carboxyl group (-COOH) at the terminus (position 1)
            and a hydroxyl (-OH) group on the opposite terminal (omega) carbon.
            Optionally, one additional hydroxyl may occur on an internal carbon.
            Any extra substituents (other than hydrogens) on the carbon chain are not allowed.
"""

from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Check whether the given SMILES corresponds to an omega-hydroxy fatty acid.
    
    Steps:
      1. Check the molecule is valid, acyclic, and composed only of C, H, and O.
      2. Confirm the presence of exactly one carboxylic acid group (SMARTS "[CX3](=O)[OX2H]").
      3. Create the subgraph of all carbon atoms. For this tree we compute its longest chain.
         (For a tree the longest path can be computed by two passes of DFS.)
      4. Confirm that the acid carbon is one of the endpoints of the longest chain.
      5. For the chain, ensure that every consecutive bond is either single or double (reject triple bonds).
      6. Check that the omega (other terminal) carbon has exactly one hydroxyl (-OH) substituent.
          Also, count any additional –OH groups on interior chain carbons; allow at most one extra.
      7. Reject if any chain carbon has any non-O, non-H substituent.
      
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): Tuple of classification and explanation.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check allowed atoms: only C (6), H (1), and O (8)
    allowed_atoms = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Molecule contains atoms other than C, H, and O"
    
    # Molecule must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected a straight-chain (acyclic) structure"

    # Work with explicit hydrogens so that we can reliably test for OH.
    mol = Chem.AddHs(mol)
    
    # 2. Identify carboxylic acid group (SMARTS: acid carbon attached to =O and -OH)
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (-COOH) group found"
    if len(acid_matches) > 1:
        return False, "Multiple carboxylic acid groups found; expected exactly one"
    # Assume the first atom in the match (the carbon) is our acid carbon.
    acid_carbon_idx = acid_matches[0][0]
    
    # 3. Build a carbon-only graph from the molecule
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    # Build dictionary: key = carbon idx, value = list of neighboring carbon indices and bond types
    carbon_graph = {}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        neighbors = []
        for bond in atom.GetBonds():
            nbr = bond.GetOtherAtom(atom)
            if nbr.GetAtomicNum() == 6:
                neighbors.append((nbr.GetIdx(), bond.GetBondType()))
        carbon_graph[idx] = neighbors
    
    # For the following, we consider the induced subgraph on carbons.
    # First, find all terminal carbons (degree 1 in the carbon graph).
    terminal_nodes = [idx for idx, nbrs in carbon_graph.items() if len(nbrs) == 1]
    if len(terminal_nodes) < 2:
        return False, "Carbon subgraph does not appear linear (cannot find enough terminal carbons)"
    
    # Because our molecule is acyclic, the carbon graph is a tree.
    # Find the longest path in the tree. We use a two-pass DFS.
    def dfs_longest(start, visited):
        # returns (farthest_node, distance, path)
        stack = [(start, 0, [start])]
        farthest_node = start
        max_dist = 0
        max_path = [start]
        while stack:
            current, dist, path = stack.pop()
            if dist > max_dist:
                max_dist = dist
                farthest_node = current
                max_path = path
            for nbr, _ in carbon_graph.get(current, []):
                if nbr not in visited:
                    visited.add(nbr)
                    stack.append((nbr, dist+1, path+[nbr]))
        return farthest_node, max_dist, max_path

    # Pick an arbitrary terminal to start
    start_node = terminal_nodes[0]
    visited = {start_node}
    node_a, _, _ = dfs_longest(start_node, visited.copy())
    # Second DFS from node_a:
    visited = {node_a}
    node_b, _, longest_path = dfs_longest(node_a, visited.copy())
    # Note: longest_path is the list of carbon indices that form the longest chain.
    
    # 4. Check that the acid carbon is an endpoint of the longest chain.
    if acid_carbon_idx not in (longest_path[0], longest_path[-1]):
        return False, "Carboxylic acid group is not located on a terminal carbon of the longest chain"
    # Designate the omega (non-acid) end:
    omega_idx = longest_path[-1] if longest_path[0] == acid_carbon_idx else longest_path[0]
    
    # 5. Check that every bond along the chain is acceptable (only single or double bonds).
    for i in range(len(longest_path)-1):
        a = longest_path[i]
        b = longest_path[i+1]
        bond = mol.GetBondBetweenAtoms(a, b)
        # Reject triple bonds (which are not expected in naturally‐occurring fatty acids)
        if bond is None:
            return False, "Chain break encountered unexpectedly"
        if bond.GetBondType() not in (Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE):
            return False, "Chain contains a bond type (e.g. triple) that is not permitted"
    
    # Helper: determine if an oxygen atom represents an -OH group.
    def is_hydroxyl(oxygen_atom):
        # oxygen should have at least one hydrogen neighbor.
        for nbr in oxygen_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                return True
        return False
    
    # 6. Now analyze substituents on the carbons in our main chain.
    chain_set = set(longest_path)
    # We now count OH substituents that are directly attached to chain carbons (but not those
    # that are part of the backbone connectivity).
    # We also require that (a) the omega carbon (the chain terminus not having the COOH group)
    # carries exactly one OH, and (b) interior chain carbons may carry at most one extra OH.
    total_chain_OH = 0
    for i, c_idx in enumerate(longest_path):
        atom = mol.GetAtomWithIdx(c_idx)
        # Get substituents that are not the chain neighbor (i.e. not in chain_set).
        off_chain = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in chain_set]
        # Also for the acid carbon we expect the -COOH group.
        if i == 0 and c_idx == acid_carbon_idx:
            # We already matched the acid -- do not count its OH (the O in COOH is part of the acid).
            # However, check that it does not have any unexpected substituents.
            for nbr in off_chain:
                # Allow oxygen only if part of COOH; we assume the acid pattern already matched.
                if nbr.GetAtomicNum() != 8:
                    return False, "Acid carbon has an unexpected substituent"
            continue
        # For the omega carbon, we require exactly one -OH substituent.
        if c_idx == omega_idx:
            oh_count = 0
            for nbr in off_chain:
                if nbr.GetAtomicNum() == 8 and is_hydroxyl(nbr):
                    oh_count += 1
                elif nbr.GetAtomicNum() != 1:  # anything that is not hydrogen is not allowed
                    return False, "Omega carbon has an unexpected substituent"
            if oh_count != 1:
                return False, "Terminal (omega) carbon does not have exactly one hydroxyl (-OH) group"
            total_chain_OH += oh_count
        else:
            # Interior carbon. They should ideally have no substituents,
            # but optionally allow one extra -OH per molecule (only on one interior carbon overall).
            oh_here = 0
            for nbr in off_chain:
                if nbr.GetAtomicNum() == 8 and is_hydroxyl(nbr):
                    oh_here += 1
                elif nbr.GetAtomicNum() != 1:
                    return False, "Interior chain carbon has an unexpected substituent (non hydroxyl)"
            total_chain_OH += oh_here

    # 7. The total count of OH groups directly attached to the chain (excluding those in the acid group)
    # must be either exactly 1 (only omega hydroxyl) or 2 (omega + one extra internal hydroxyl)
    if total_chain_OH not in (1, 2):
        return False, f"Number of hydroxyl substituents attached to the chain is {total_chain_OH}; expected 1 or 2"
        
    return True, "Molecule is a straight-chain fatty acid with a COOH at a terminal carbon and an -OH at the omega position (with at most one extra internal hydroxyl allowed)"

    
# Example test (you may remove or modify as needed)
if __name__ == '__main__':
    # Test with one example: (2E)-18-hydroxyoctadec-2-enoic acid
    test_smiles = "OCCCCCCCCCCCCCCC\\C=C\\C(O)=O"
    result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
    print(result, reason)