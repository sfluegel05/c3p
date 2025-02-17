"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol

An alditol is an acyclic polyol having the general formula:
   HOCH2[CH(OH)]nCH2OH 
(formally derivable from an aldose by reduction of its carbonyl group).

In the previous approach the identification of candidate carbons depended on a helper function 
(checking for an alcohol) that was overly strict. The improved version redefines is_alcohol to:
  • Check that the atom is an oxygen.
  • Verify that every bond from the oxygen is a single bond.
  • Check that the oxygen is attached to exactly one non‐hydrogen atom.
  • Ensure that the oxygen is attached to at least one hydrogen.
Candidate carbons are those sp³ carbons (atomic number 6) that are not in a ring and have at least 
one attached free hydroxyl. We then build a connectivity graph among candidates, linking carbons connected 
by single bonds (outside rings). Finally, we search for any contiguous chain (via DFS) of at least three 
candidate carbons. If such a chain is found, the molecule is classified as an alditol.
"""
from rdkit import Chem
from collections import deque

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    
    An alditol is defined as an acyclic polyol conforming to the formula HOCH2[CH(OH)]nCH2OH.
    We identify candidate carbon atoms (sp³ and not in a ring) that bear a free hydroxyl group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        (bool, str): Tuple where the boolean is True if an alditol chain is found,
                     and False otherwise. The string provides the reasoning.
    """
    # Parse SMILES string and add explicit hydrogens
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Improved helper function: determine if a given oxygen atom is a free hydroxyl.
    def is_alcohol(atom):
        # Must be oxygen.
        if atom.GetSymbol() != "O":
            return False
        # Verify that all bonds from oxygen are single.
        for bond in atom.GetBonds():
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                return False
        # Check that it is attached to exactly one heavy (non‐hydrogen) atom.
        heavy_neighbors = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            return False
        # Must be attached to at least one hydrogen (explicitly present).
        hydrogen_neighbors = [nb for nb in atom.GetNeighbors() if nb.GetAtomicNum() == 1]
        if len(hydrogen_neighbors) < 1:
            return False
        return True

    candidate_indices = []
    # Identify candidate carbon atoms based on:
    #  • Carbon (atomic number 6)
    #  • sp³ hybridized (to favor acyclic centers)
    #  • Not in a ring (acyclic chain is desired)
    #  • Has at least one attached free hydroxyl group.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        if atom.IsInRing():
            continue
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
    # Two candidate carbons are connected if they share a single bond that is not part of a ring.
    candidate_set = set(candidate_indices)
    graph = {idx: [] for idx in candidate_set}
    for idx in candidate_indices:
        atom = mol.GetAtomWithIdx(idx)
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in candidate_set:
                bond = mol.GetBondBetweenAtoms(idx, nb_idx)
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE and not bond.IsInRing():
                    if nb_idx not in graph[idx]:
                        graph[idx].append(nb_idx)
                    if idx not in graph[nb_idx]:
                        graph[nb_idx].append(idx)
    
    # First, get connected components from the candidate graph.
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
    
    # For each connected component, use DFS to search for a chain of at least three candidate carbons.
    def dfs_paths(g, current, visited):
        yield visited
        for nbr in g[current]:
            if nbr not in visited:
                yield from dfs_paths(g, nbr, visited + [nbr])
    
    for comp in components:
        if len(comp) < 3:
            continue
        # Build subgraph for the current component.
        comp_graph = {node: graph[node] for node in comp}
        for start in comp_graph:
            for path in dfs_paths(comp_graph, start, [start]):
                if len(path) >= 3:
                    return True, f"Found alditol chain of {len(path)} contiguous candidate carbons bearing free –OH groups."
    
    return False, "No linear alditol chain pattern detected"

# Example usage and testing
if __name__ == "__main__":
    # A selection of SMILES examples (both positive and negatives) are provided for testing purposes.
    test_smiles = [
        "C([C@H]([C@H]([C@@H]([C@H]([13CH2]O)O)O)O)O)O",  # D-(1-(13)C)glucitol
        "OC[C@@H](O)[C@H](O)CO",                           # D-threitol
        "C[C@@](O)(CO)[C@@H](O)CO",                        # 2-methylerythritol
        "OC[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO",         # L-iditol
        "OC[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO",         # D-iditol
        "OC[C@H](O)[C@H](O)[C@H](O)[C@@H](O)CO",          # L-altritol
        "OCC(O)CO"                                       # glycerol (not an alditol chain in our context)
    ]
    for s in test_smiles:
        res, msg = is_alditol(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n{'-'*60}")