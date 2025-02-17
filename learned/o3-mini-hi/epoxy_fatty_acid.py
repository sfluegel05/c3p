"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
Improved criteria:
  - Contains a terminal carboxylic acid group (C(=O)[O;H]) where the acid carbon is attached 
    to exactly one other carbon.
  - Has a long acyclic chain (i.e. a chain of non‐ring carbon atoms) originating from the carbon 
    attached to the terminal acid; we require that the longest such path contains at least 12 carbons.
  - The structure should have a high ratio of acyclic carbons relative to all carbons (>0.5).
  - Contains exactly one epoxide ring defined as a three-membered ring with two carbons and one oxygen (SMARTS: [CX2]1[OX2][CX2]1).
  - The overall molecular weight must exceed 200 Da.
  
Note: Extra rings other than the epoxide are allowed provided that the longest acyclic chain is long enough.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from collections import deque

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    
    The criteria are:
      - It must have a terminal carboxylic acid group (C(=O)[O;H]) where the carbonyl carbon is attached to only one carbon.
      - The neighbor carbon bonded to the acid group must connect to a sufficiently long chain (at least 12 acyclic carbons in a connected chain).
      - The molecule should have a high ratio of acyclic (non-ring) carbons (>0.5 of all carbon atoms).
      - It must contain exactly one epoxide ring defined by the SMARTS "[CX2]1[OX2][CX2]1".
      - Its molecular weight is above 200 Da.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the epoxy fatty acid criteria; False otherwise.
        str: A message stating the reason for the classification decision.
    """
    
    # Parse SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Identify terminal carboxylic acid group.
    # The SMARTS "C(=O)[O;H]" matches a carboxyl group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found to signify a fatty acid"
    
    terminal_acid_found = False
    acid_neighbor_idx = None
    # Loop over matches to find one where the acid carbon (index at position 0) has exactly one carbon neighbor.
    for match in acid_matches:
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Consider only neighbors that are carbon atoms.
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            # Found a terminal acid: the acid carbon is attached to exactly one carbon.
            terminal_acid_found = True
            acid_neighbor_idx = carbon_neighbors[0]
            break
    
    if not terminal_acid_found:
        return False, "Carboxylic acid group found, but it is not terminal (acid carbon attached to >1 carbon)"
    
    # Check basic fatty acid properties: overall number of carbons and molecular weight.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if total_carbons < 12:
        return False, f"Too few carbon atoms ({total_carbons}) to be a fatty acid"
    
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.2f} Da) too low for a fatty acid"
    
    # Check acyclic (non-ring) carbon ratio.
    acyclic_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing())
    acyclic_ratio = acyclic_carbons / total_carbons if total_carbons > 0 else 0
    if acyclic_ratio < 0.5:
        return False, f"Low acyclic carbon ratio ({acyclic_ratio:.2f}); structure too cyclic to be a fatty acid"
    
    # From the acid neighbor, build a graph of acyclic carbons.
    # Only include carbon atoms that are not in any ring.
    acyclic_atoms = set(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing())
    
    if acid_neighbor_idx not in acyclic_atoms:
        return False, "Acid neighbor is not part of an acyclic chain"
    
    # Build undirected graph from the acyclic carbon atoms.
    graph = {idx: [] for idx in acyclic_atoms}
    for idx in acyclic_atoms:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() in acyclic_atoms):
                graph[idx].append(nbr.GetIdx())
    
    # Use BFS to determine the longest distance from the starting atom (acid neighbor) in this acyclic graph.
    distances = {acid_neighbor_idx: 0}
    queue = deque([acid_neighbor_idx])
    while queue:
        current = queue.popleft()
        for neighbor in graph[current]:
            if neighbor not in distances:
                distances[neighbor] = distances[current] + 1
                queue.append(neighbor)
    
    # The chain length is measured as the maximum number of carbons from the acid neighbor along the acyclic graph.
    # We add 1 to count the starting acid neighbor carbon.
    if distances:
        chain_length = max(distances.values()) + 1
    else:
        chain_length = 1  # Only the acid neighbor
    
    if chain_length < 12:
        return False, f"Longest acyclic carbon chain from the acid group is too short (length: {chain_length})"
    
    # Identify epoxide rings: a 3-membered ring with one oxygen and two carbons.
    epoxide_pattern = Chem.MolFromSmarts("[CX2]1[OX2][CX2]1")
    epoxide_matches = mol.GetSubstructMatches(epoxide_pattern)
    # Use a set of sorted tuples to count unique epoxide rings.
    epoxide_sets = set(tuple(sorted(match)) for match in epoxide_matches)
    if len(epoxide_sets) != 1:
        return False, f"Expected exactly one epoxide ring, found {len(epoxide_sets)}"
    
    return True, "Contains terminal carboxylic acid group, a long acyclic chain, and one epoxide ring indicative of an epoxy fatty acid"

# Example usage (uncomment to test):
# test_smiles = "C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O"  # (5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid
# result, reason = is_epoxy_fatty_acid(test_smiles)
# print(result, reason)