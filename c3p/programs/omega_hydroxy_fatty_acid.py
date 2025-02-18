"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
"""
Classifies: Omega-hydroxy fatty acid
Definition: A naturally-occurring straight-chain fatty acid (acyclic and unbranched)
            that contains only C, H, and O, with a terminal carboxyl group (-COOH)
            at position 1 and a hydroxyl (-OH) group on the opposite terminal (omega) carbon.
"""

from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if the given molecule (as a SMILES string) is an omega-hydroxy fatty acid.
    The criteria are:
      1. The molecule contains only C, H, and O.
      2. The molecule is acyclic.
      3. There is exactly one carboxylic acid (-COOH) group (identified via SMARTS "[CX3](=O)[OX2H]").
      4. When considering only its carbon atoms the molecule forms a single unbranched (linear) chain:
         exactly two carbon atoms have only one neighboring carbon and all other carbons have two.
      5. One chain terminus must be the acid carbon and the other (omega) carbon must bear an -OH group.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the omega-hydroxy fatty acid criteria; False otherwise.
        str: Explanation of the result.
    """
    # Parse the molecule from SMILES.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check that only C, H, and O (atomic numbers: 6, 1, 8) are present.
    allowed_atoms = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atoms:
            return False, "Molecule contains atoms other than C, H, and O"
    
    # 2. Exclude cyclic molecules; the fatty acid must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic; expected a straight-chain (acyclic) structure"
    
    # It is helpful for recognising -OH groups to include hydrogens explicitly.
    mol = Chem.AddHs(mol)
    
    # 3. Identify the carboxylic acid group using SMARTS.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid (-COOH) group found"
    # For a fatty acid we require exactly one -COOH group.
    if len(acid_matches) > 1:
        return False, "Multiple carboxylic acid groups found; fatty acid should contain just one"
    # We assume the first atom in the match is the acid carbon.
    acid_carbon_idx = acid_matches[0][0]
    
    # 4. Build an induced graph using only carbon atoms.
    # Get all carbon atoms (by their indices) in the molecule.
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    
    # Create a dictionary mapping each carbon index to its neighboring carbon indices.
    carbon_graph = {}
    for idx in carbon_indices:
        atom = mol.GetAtomWithIdx(idx)
        # Only count neighbors that are carbons.
        neighbors = [nbr.GetIdx() for nbr in atom.GetNeighbors() if nbr.GetAtomicNum() == 6]
        carbon_graph[idx] = neighbors
    
    # Check that the induced carbon subgraph is connected.
    # We'll start from an arbitrary carbon (say, the acid carbon) and do a DFS.
    visited = set()
    stack = [acid_carbon_idx]
    while stack:
        current = stack.pop()
        if current not in visited:
            visited.add(current)
            stack.extend(carbon_graph[current])
    if set(carbon_indices) != visited:
        return False, "Not all carbon atoms are connected in a single chain"
    
    # For a perfect linear chain:
    # - exactly two carbons (the terminals) should have degree 1;
    # - all others must have degree 2.
    terminal_nodes = [idx for idx, nbrs in carbon_graph.items() if len(nbrs) == 1]
    for idx, nbrs in carbon_graph.items():
        if idx not in terminal_nodes and len(nbrs) != 2:
            return False, "Carbon subgraph is branched; expected a straight-chain fatty acid"
    if len(terminal_nodes) != 2:
        return False, "Carbon subgraph does not have exactly two termini; found %d" % len(terminal_nodes)
    
    # 5. One terminal must be the acid carbon.
    if acid_carbon_idx not in terminal_nodes:
        return False, "Carboxylic acid group is not located on a terminal carbon"
    # Identify the omega carbon as the other terminal.
    omega_candidate_idx = terminal_nodes[0] if terminal_nodes[1] == acid_carbon_idx else terminal_nodes[1]
    omega_candidate = mol.GetAtomWithIdx(omega_candidate_idx)
    
    # Check that the omega candidate bears a hydroxyl (-OH) substituent.
    has_oh = False
    for nbr in omega_candidate.GetNeighbors():
        # Only consider oxygen atoms.
        if nbr.GetAtomicNum() == 8:
            # We require the oxygen to have at least one hydrogen attached.
            # (This distinguishes a hydroxyl from a carbonyl oxygen.)
            if any(neigh.GetAtomicNum() == 1 for neigh in nbr.GetNeighbors()):
                has_oh = True
                break
    if not has_oh:
        return False, "Terminal (omega) carbon lacks a hydroxyl (-OH) group"
    
    return True, "Molecule is a straight-chain fatty acid with a COOH at C1 and an -OH at the omega position"

    
# Example test call (you can modify or remove this as needed):
if __name__ == '__main__':
    test_smiles = "OCCCCCCCCCCCCCCC\\C=C\\C(O)=O"  # (2E)-18-hydroxyoctadec-2-enoic acid
    result, reason = is_omega_hydroxy_fatty_acid(test_smiles)
    print(result, reason)