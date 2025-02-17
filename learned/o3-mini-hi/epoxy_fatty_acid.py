"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
Criteria:
  - Must have a terminal carboxylic acid group (C(=O)[O;H]), wherein the acid carbon is attached to exactly one carbon.
  - From the acid, the connected acyclic (non‐ring) carbon chain must be long (at least 12 carbons).
  - The ratio of acyclic carbons (outside rings) to total carbons must exceed 0.5.
  - Contains exactly one epoxide ring, defined as a three‐membered ring that has exactly 2 carbons and 1 oxygen.
  - The overall molecular weight must be above 200 Da.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    
    The criteria are:
      - Contains a terminal carboxylic acid group (matched by C(=O)[O;H]) where the acid carbon
        is attached to exactly one carbon.
      - The neighbor carbon (attached to the acid group) must lead to an acyclic chain of at
        least 12 carbons, computed via a DFS over all acyclic (non‐ring) carbon atoms.
      - The ratio of acyclic carbons (non‐ring carbons) to total carbons must exceed 0.5.
      - Contains exactly one epoxide ring. For robustness, we count every ring in the molecule
        which is of size 3 having exactly two carbons and one oxygen.
      - The overall molecular weight is above 200 Da.
    
    Returns:
        (bool): True if the molecule meets the epoxy fatty acid criteria; False otherwise.
        (str): Reason for classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # 1. Check for the terminal carboxylic acid group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found to signify a fatty acid"
    
    terminal_acid_found = False
    acid_neighbor_idx = None
    # Look for an acid match where the acid carbon is attached to exactly one carbon.
    for match in acid_matches:
        acid_carbon = mol.GetAtomWithIdx(match[0])
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            terminal_acid_found = True
            acid_neighbor_idx = carbon_neighbors[0]
            break
    if not terminal_acid_found:
        return False, "Carboxylic acid group found, but it is not terminal (acid carbon attached to >1 carbon)"
    
    # 2. Check basic fatty acid properties.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 12:
        return False, f"Too few carbon atoms ({total_carbons}) to be a fatty acid"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.2f} Da) too low for a fatty acid"
    
    # 3. Check the ratio of acyclic carbons.
    acyclic_carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if total_carbons == 0:
        return False, "No carbon atoms found"
    acyclic_ratio = len(acyclic_carbons) / total_carbons
    if acyclic_ratio < 0.5:
        return False, f"Low acyclic carbon ratio ({acyclic_ratio:.2f}); structure too cyclic to be a fatty acid"
    
    # 4. Determine the longest acyclic chain starting from the acid neighbor.
    # Build a graph of acyclic (non‐ring) carbon atoms. The nodes are the atom indices.
    acyclic_atoms = set(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing())
    if acid_neighbor_idx not in acyclic_atoms:
        return False, "Acid neighbor is not part of an acyclic chain"
    
    # Build an adjacency list graph.
    graph = {idx: [] for idx in acyclic_atoms}
    for idx in acyclic_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() in acyclic_atoms:
                graph[idx].append(nbr.GetIdx())
                
    # Now use DFS to compute the longest simple path starting from acid_neighbor_idx.
    def dfs(node, visited):
        max_length = 1  # Count current node.
        for neighbor in graph[node]:
            if neighbor not in visited:
                length = 1 + dfs(neighbor, visited | {neighbor})
                if length > max_length:
                    max_length = length
        return max_length
    
    chain_length = dfs(acid_neighbor_idx, {acid_neighbor_idx})
    if chain_length < 12:
        return False, f"Longest acyclic carbon chain from the acid group is too short (length: {chain_length})"
    
    # 5. Count epoxide rings.
    # Instead of relying solely on a SMARTS, we iterate over all rings of size 3 and count those with exactly 2 carbons and 1 oxygen.
    ring_info = mol.GetRingInfo()
    epoxide_count = 0
    for ring in ring_info.AtomRings():
        if len(ring) == 3:
            atom_nums = [mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring]
            # Check if sorted atomic numbers match two carbons (6) and one oxygen (8).
            if sorted(atom_nums) == [6, 6, 8]:
                epoxide_count += 1
    if epoxide_count != 1:
        return False, f"Expected exactly one epoxide ring, found {epoxide_count}"
    
    return True, ("Contains terminal carboxylic acid group with a sufficiently long acyclic chain (length: "
                  f"{chain_length}), a high acyclic carbon ratio, and one epoxide ring indicative of an epoxy fatty acid")

# Example usage (uncomment to test):
# test_smiles = "C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O"  # (5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid
# result, reason = is_epoxy_fatty_acid(test_smiles)
# print(result, reason)