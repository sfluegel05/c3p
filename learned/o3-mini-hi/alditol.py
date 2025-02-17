"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol

An alditol is an acyclic polyol having the general formula:
  HOCH2[CH(OH)]nCH2OH
(i.e. formally derivable from an aldose by reduction of its carbonyl group).

This implementation first adds explicit hydrogens, then identifies candidate carbon
atoms. A candidate carbon is defined as an sp3 (saturated) carbon not in a ring and that
bears at least one “free” –OH group. Here, a free –OH is defined as an oxygen atom that:
  • has atomic number 8,
  • is attached to exactly one heavy atom (thus not bridging two carbons), and
  • carries at least one hydrogen.
After collecting candidate carbons (presumed to be part of the main open-chain portion)
a graph is built by connecting candidate carbons that share a single, non‐ring bond.
Finally, we search for a simple (acyclic) contiguous chain of at least three candidate carbons.
If found, the molecule is classified as an alditol.
"""

from rdkit import Chem
from collections import deque

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): A tuple where the first element is True if an alditol chain is found,
                   and False otherwise. The second element gives the reasoning.
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to help account for –OH groups
    mol = Chem.AddHs(mol)
    
    # Define a helper function: determines if an oxygen atom is a free hydroxyl.
    # A free hydroxyl (unmodified –OH) is defined here as:
    #   • An oxygen (atomic number 8)
    #   • That is attached to exactly one heavy atom (i.e. not bridging two carbons)
    #   • None of its bonds are double bonds (so that it is not a carbonyl oxygen)
    #   • And it carries at least one hydrogen.
    def is_alcohol(oxygen_atom):
        if oxygen_atom.GetAtomicNum() != 8:
            return False
        # Count heavy (non-hydrogen) neighbors.
        heavy_neighbors = sum(1 for nb in oxygen_atom.GetNeighbors() if nb.GetAtomicNum() != 1)
        if heavy_neighbors != 1:
            return False
        # Verify that no bond to this oxygen is a double bond.
        for bond in oxygen_atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                return False
        # Check that the oxygen carries at least one hydrogen.
        if oxygen_atom.GetTotalNumHs() < 1:
            return False
        return True

    # Identify candidate carbons. We require:
    #   • The atom is carbon (atomic number 6)
    #   • It is sp3 hybridized
    #   • It is not in a ring (acyclic chain is desired)
    #   • It has a neighboring oxygen that qualifies as a free –OH.
    candidate_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        if atom.IsInRing():
            continue
        # Check if at least one neighbor is a free hydroxyl oxygen.
        has_free_oh = False
        for nb in atom.GetNeighbors():
            if is_alcohol(nb):
                has_free_oh = True
                break
        if has_free_oh:
            candidate_indices.append(atom.GetIdx())
    
    if not candidate_indices:
        return False, "No candidate carbon bearing a free –OH group found"
    
    # Build a connectivity graph among candidate carbons.
    # We define an edge if two candidate carbons are directly bonded by a single bond
    # that is not part of a ring.
    candidate_set = set(candidate_indices)
    graph = {idx: [] for idx in candidate_set}
    for idx in candidate_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in candidate_set:
                bond = mol.GetBondBetweenAtoms(idx, nb_idx)
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    if not bond.IsInRing():
                        if nb_idx not in graph[idx]:
                            graph[idx].append(nb_idx)
                        if idx not in graph[nb_idx]:
                            graph[nb_idx].append(idx)
    
    # Get connected components of the candidate graph.
    def get_connected_components(g):
        seen = set()
        components = []
        for node in g:
            if node not in seen:
                comp = set()
                queue = deque([node])
                seen.add(node)
                while queue:
                    curr = queue.popleft()
                    comp.add(curr)
                    for nbr in g[curr]:
                        if nbr not in seen:
                            seen.add(nbr)
                            queue.append(nbr)
                components.append(comp)
        return components
    
    components = get_connected_components(graph)
    
    # In each connected component, try to find a simple acyclic path 
    # (via depth-first search) that contains at least 3 candidate carbons.
    def dfs_paths(g, current, visited):
        yield visited
        for nbr in g[current]:
            if nbr not in visited:
                yield from dfs_paths(g, nbr, visited + [nbr])
    
    for comp in components:
        if len(comp) < 3:
            continue
        comp_graph = {node: graph[node] for node in comp}
        for start in comp_graph:
            for path in dfs_paths(comp_graph, start, [start]):
                if len(path) >= 3:
                    return True, f"Found alditol chain of {len(path)} contiguous candidate carbons bearing free –OH groups."
    
    return False, "No linear alditol chain pattern detected"

# Example usage and testing – you can run this module to check several examples:
if __name__ == "__main__":
    test_smiles = [
        "C([C@H]([C@H]([C@@H]([C@H]([13CH2]O)O)O)O)O)O",  # D-(1-(13)C)glucitol
        "OC[C@@H](O)[C@H](O)CO",                           # D-threitol
        "C[C@@](O)(CO)[C@@H](O)CO",                        # 2-methylerythritol
        "OC[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO",         # L-iditol
        "O([C@@H]([C@H](O[C@@H]1O[C@H]([C@@H](O)[C@@H](O)[C@@H]1O)C)[C@@H](O)CO)[C@H](O)CO)[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO",  # complex glycoside example
        "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)C1OC(COC2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)CO)CO",  # another glycoside
        "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OCC(O)CO)CO",  # 1-O-(beta-D-glucosyl)glycerol variant
        "C(CO)(CO)O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO",  # 2-O-(beta-D-glucosyl)glycerol variant
        "OC[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO",         # D-iditol
        "OC[C@H](O)[C@H](O)[C@H](O)[C@@H](O)CO",          # L-altritol
        "OCC(O)CO"                                       # glycerol, smallest alditol
    ]
    for s in test_smiles:
        res, msg = is_alditol(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n{'-'*60}")