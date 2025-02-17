"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: Omega-hydroxy fatty acid
Definition: A naturally-occurring straight-chain fatty acid (acyclic, unbranched)
            composed only of C, H, O that contains one carboxyl group (-COOH) at C1
            and a hydroxyl (-OH) group on the opposite terminal (omega) carbon.
            Optionally, one additional hydroxyl may occur on an internal carbon.
            (Sugar acids and cyclic molecules are rejected.)
"""

from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if the given molecule (as a SMILES string) is an omega-hydroxy fatty acid.
    
    The criteria are:
      1. The molecule contains only C, H, and O.
      2. The molecule is acyclic.
      3. There is exactly one carboxylic acid (-COOH) group (identified via SMARTS "[CX3](=O)[OX2H]").
      4. When considering the carbon atoms, they form a single unbranched chain.
         (That is, exactly two carbons (the termini) have only one neighbor within the set.)
      5. The acid carbon (of the -COOH group) must be one terminal carbon.
      6. The non-acid terminal ("omega carbon") must bear exactly one hydroxyl (-OH) substituent.
      7. Aside from the omega OH, at most one interior chain carbon may bear an extra OH.
         (This helps to reject polyhydroxylated sugar acids.)
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the omega-hydroxy fatty acid criteria; False otherwise.
        str: Explanation of the result.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check that only C, H, and O are present.
    allowed_atoms = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Molecule contains atoms other than C, H, and O"
    
    # 2. Molecule must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected a straight-chain (acyclic) structure"
    
    # It will be useful to have explicit hydrogens.
    mol = Chem.AddHs(mol)
    
    # 3. Identify the carboxylic acid group.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (-COOH) group found"
    if len(acid_matches) > 1:
        return False, "Multiple carboxylic acid groups found; fatty acid should contain just one"
    # We assume the first atom in the match is the acid carbon.
    acid_carbon_idx = acid_matches[0][0]
    
    # 4. Build the induced carbon graph.
    # Get indices for all carbon atoms.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    carbon_graph = {}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        # Count only carbon neighbors.
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        carbon_graph[idx] = neighbors
    
    # Check connectivity in the carbon subgraph.
    visited = set()
    stack = [acid_carbon_idx]
    while stack:
        current = stack.pop()
        if current not in visited:
            visited.add(current)
            stack.extend(carbon_graph[current])
    if set(carbon_indices) != visited:
        return False, "Not all carbon atoms are connected in a single chain"
    
    # In a linear (unbranched) chain, exactly two carbons have degree 1.
    terminal_nodes = [idx for idx, nbrs in carbon_graph.items() if len(nbrs) == 1]
    for idx, nbrs in carbon_graph.items():
        if idx not in terminal_nodes and len(nbrs) != 2:
            return False, "Carbon subgraph is branched; expected a straight-chain fatty acid"
    if len(terminal_nodes) != 2:
        return False, "Carbon subgraph does not have exactly two termini; found %d" % len(terminal_nodes)
    
    # 5. The acid carbon must be one terminal.
    if acid_carbon_idx not in terminal_nodes:
        return False, "Carboxylic acid group is not located on a terminal carbon"
    # Designate the other terminal as the omega candidate.
    omega_idx = terminal_nodes[0] if terminal_nodes[1] == acid_carbon_idx else terminal_nodes[1]
    
    # 6. Order the carbon chain starting from the acid carbon.
    chain = [acid_carbon_idx]
    current = acid_carbon_idx
    prev = -1
    # Because the chain is unbranched there is exactly one neighbor (other than prev) each time.
    while current != omega_idx:
        # Get the neighbor in the chain that is not the previous carbon.
        nxt = [nbr for nbr in carbon_graph[current] if nbr != prev]
        if len(nxt) != 1:
            return False, "Unable to determine unique linear ordering of the carbon chain"
        nxt = nxt[0]
        chain.append(nxt)
        prev, current = current, nxt
    # Now chain[0] is acid carbon and chain[-1] is omega carbon.
    
    # 7. Check substituents on chain carbons.
    # Define a helper: a neighbor is considered an -OH substituent if it is an oxygen
    # that in turn has at least one hydrogen.
    def is_hydroxyl(os_atom):
        # os_atom is expected to be an oxygen.
        for nbr in os_atom.GetNeighbors():
            if nbr.GetAtomicNum() == 1:
                return True
        return False

    interior_OH_count = 0  # count OH's on carbons excluding acid and omega.
    
    # For each carbon in the chain:
    for i, c_idx in enumerate(chain):
        atom = mol.GetAtomWithIdx(c_idx)
        # Get all neighbors that are not part of the carbon chain:
        non_chain_neighbors = [nbr for nbr in atom.GetNeighbors() if nbr.GetIdx() not in chain and nbr.GetAtomicNum() != 1]
        
        # For the acid carbon (chain[0]), we expect the -COOH group.
        if i == 0:
            # The acid carbon is already confirmed by SMARTS.
            continue
        # For the omega carbon, check that it has exactly one -OH substituent.
        if i == len(chain) - 1:
            oh_count = 0
            for nbr in non_chain_neighbors:
                if nbr.GetAtomicNum() == 8 and is_hydroxyl(nbr):
                    oh_count += 1
            if oh_count != 1:
                return False, "Terminal (omega) carbon does not have exactly one hydroxyl (-OH) group"
            continue
        # For an interior carbon, ideally there are no extra substituents.
        # However, some fatty acids may have one extra -OH.
        oh_here = 0
        for nbr in non_chain_neighbors:
            # If attached atom is oxygen and displays a hydroxyl, count it.
            if nbr.GetAtomicNum() == 8 and is_hydroxyl(nbr):
                oh_here += 1
            else:
                # If there is any non-oxygen substituent (or oxygen not behaving like -OH),
                # then the chain is not a simple alkyl backbone.
                return False, "Interior chain carbon has an unexpected substituent (non hydroxyl)"
        if oh_here > 1:
            return False, "An interior carbon has >1 hydroxyl substituents (unexpected)"
        interior_OH_count += oh_here

    # Allow at most one interior hydroxyl (so total OH=omega + interior is either 1 or 2).
    if interior_OH_count > 1:
        return False, f"Too many hydroxyl substituents on the chain; found {1 + interior_OH_count} total"
    
    return True, "Molecule is a straight-chain fatty acid with a COOH at C1 and an -OH at the omega position"

    
# Example test (you may remove or modify as needed)
if __name__ == '__main__':
    test_smiles = "OCCCCCCCCCCCCCCC\\C=C\\C(O)=O"  # (2E)-18-hydroxyoctadec-2-enoic acid
    result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
    print(result, reason)