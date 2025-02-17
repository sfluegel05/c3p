"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: Epoxy fatty acid
Definition: A heterocyclic fatty acid containing an epoxide ring as part of its structure.
Criteria:
  - Must have a terminal carboxylic acid group (C(=O)[O;H]) in which the acid carbon is attached to exactly one carbon.
  - From that acid group, the connected acyclic (non‐ring) carbon component must be long (at least 12 carbons).
  - The ratio of acyclic (non‐ring) carbons to total carbons (in the full molecule) must exceed 0.5.
  - Contains exactly one epoxide ring. We define an epoxide ring as a 3‐membered ring (in the molecule’s ring info)
    having exactly two carbons and one oxygen – and which is part of the fatty acyl chain.
  - The overall molecular weight must be above 200 Da.
  
If any of these criteria fail (or if the SMILES cannot be parsed), the function returns False with a reason.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.
    
    Steps/Criteria:
      1. Parse the molecule and verify its molecular weight > 200 Da.
      2. Find a terminal carboxylic acid group (match: C(=O)[O;H]) in which the acid carbon is attached to exactly one carbon.
      3. Build the graph of acyclic (non‐ring) carbon atoms and, starting from the acid neighbor,
         compute:
           (a) the largest connected acyclic carbon subcomponent (via BFS) – so that we later force the epoxide ring to be in that chain.
           (b) the longest simple path (via DFS) from that acid neighbor within the acyclic subgraph.
         One candidate acid group is acceptable if its longest acyclic chain length is at least 12 carbons.
      4. Verify that overall the ratio of acyclic carbons (in the molecule) to total carbons exceeds 0.5.
      5. Count epoxide rings – only counting those rings (from mol.GetRingInfo()) of size 3 that consist of exactly two carbons and one oxygen, and that are fully contained in the fatty acyl chain (the connected acyclic subcomponent).
         Exactly one such epoxide ring is required.
    
    Returns:
        (bool): True if all criteria are met, else False.
        (str): Reason for the classification.
    """
    # Parse the molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check overall molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, f"Molecular weight ({mol_wt:.2f} Da) too low for a fatty acid"
    
    # Find all carboxylic acid groups that match our pattern.
    acid_pattern = Chem.MolFromSmarts("C(=O)[O;H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found to signify a fatty acid"
    
    # We will try all acid matches; acceptance of one candidate is enough.
    candidate_found = False
    candidate_chain_length = 0
    candidate_reason = ""
    
    # Build overall count of carbon atoms.
    total_carbons = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    # Also count acyclic carbons (non-ring carbons).
    acyclic_atoms_global = set(atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing())
    if total_carbons == 0:
        return False, "No carbon atoms found in the molecule"
    acyclic_ratio = len(acyclic_atoms_global) / total_carbons
    if acyclic_ratio < 0.5:
        return False, f"Low acyclic carbon ratio ({acyclic_ratio:.2f}); structure too cyclic to be a fatty acid"
    
    # Build graph for acyclic carbon atoms.
    # Nodes: indices of carbon atoms not in rings.
    acyclic_graph = {idx: [] for idx in acyclic_atoms_global}
    for idx in acyclic_atoms_global:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            if nbr.GetAtomicNum() == 6 and (nbr.GetIdx() in acyclic_atoms_global):
                acyclic_graph[idx].append(nbr.GetIdx())
    
    # DFS to compute the longest simple path from a given starting node in the acyclic graph.
    def dfs_longest_path(node, visited):
        max_length = 1  # count the current node
        for neighbor in acyclic_graph.get(node, []):
            if neighbor not in visited:
                path_length = 1 + dfs_longest_path(neighbor, visited | {neighbor})
                if path_length > max_length:
                    max_length = path_length
        return max_length

    # For epoxide search:
    # We'll later get the connected set (component) of acyclic carbons attached to the acid neighbor.
    def bfs_component(start):
        comp = {start}
        queue = [start]
        while queue:
            current = queue.pop(0)
            for nbr in acyclic_graph.get(current, []):
                if nbr not in comp:
                    comp.add(nbr)
                    queue.append(nbr)
        return comp

    # Loop over the acid matches to find one acceptable candidate.
    candidate_chain_lengths = []
    for match in acid_matches:
        acid_carbon = mol.GetAtomWithIdx(match[0])  # acid carbon (the C in C(=O)[O;H])
        # Collect carbon neighbors of the acid carbon.
        carbon_neighbors = [nbr.GetIdx() for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) != 1:
            continue  # Not terminal.
        acid_neighbor = carbon_neighbors[0]
        # Only continue if the acid_neighbor is acyclic.
        if acid_neighbor not in acyclic_atoms_global:
            continue
        
        # Compute the longest simple chain (number of carbons) starting from acid_neighbor.
        chain_length = dfs_longest_path(acid_neighbor, {acid_neighbor})
        # Also get the connected acyclic component from this neighbor.
        chain_component = bfs_component(acid_neighbor)
        
        candidate_chain_lengths.append(chain_length)
        if chain_length >= 12:
            candidate_found = True
            candidate_chain_length = chain_length
            # Now check epoxide rings that lie in the fatty acyl chain.
            ring_info = mol.GetRingInfo()
            epoxide_count = 0
            for ring in ring_info.AtomRings():
                if len(ring) != 3:
                    continue
                # Get the atomic numbers of atoms in this ring.
                atom_nums = sorted(mol.GetAtomWithIdx(idx).GetAtomicNum() for idx in ring)
                if atom_nums != [6, 6, 8]:
                    continue
                # Now, require that the entire ring is part of our acyclic fatty acid chain.
                # (This ensures the epoxide is “in” the fatty acyl chain.)
                if set(ring).issubset(chain_component):
                    epoxide_count += 1
            if epoxide_count != 1:
                candidate_reason = f"Expected exactly one epoxide ring within the fatty acyl chain, found {epoxide_count}"
                continue  # try next candidate if any
            # If we reach here, all criteria (so far) are met.
            candidate_reason = ("Contains terminal carboxylic acid group with a sufficiently long acyclic chain (length: "
                                f"{chain_length}), a high acyclic carbon ratio ({acyclic_ratio:.2f}), and one epoxide ring indicative of an epoxy fatty acid")
            break

    if not candidate_found:
        # If we had candidates but none with chain length >= 12 then report the max found.
        if candidate_chain_lengths:
            max_chain = max(candidate_chain_lengths)
            return False, f"Longest acyclic carbon chain from the acid group is too short (length: {max_chain})"
        else:
            return False, "No terminal carboxylic acid group with a proper acyclic attachment found"
    
    return True, candidate_reason

# Example usage:
# test_smiles = "C(CCC/C=C\\C[C@@H]1/C(/O1)=C/C=C\\C/C=C\\CCCCC)(=O)O"  # (5Z,8R,9Z,11Z,14Z)-8,9-epoxyicosatetraenoic acid
# result, reason = is_epoxy_fatty_acid(test_smiles)
# print(result, reason)