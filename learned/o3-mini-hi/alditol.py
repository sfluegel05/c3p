"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol

An alditol is an acyclic polyol having the general formula:
  HOCH2[CH(OH)]nCH2OH
(i.e. it is formally derivable from an aldose by reduction of its carbonyl group).

This implementation first adds explicit hydrogens to the molecule.
Then it identifies candidate carbon atoms that are sp³, not in a ring, and 
that bear at least one “free” (alcohol) –OH group. (Here, “free –OH” is defined 
as an oxygen atom that is singly bound (no double‐bond) and that carries at least 
one hydrogen.)
Next, a graph is built from these candidate carbons (an edge is drawn if two 
candidate carbons share a single, non‐ring bond). Finally, we scan each connected 
component for a simple linear (acyclic) path of at least three continuous candidate 
carbons. If such a chain is found, we return True with details, else False.
"""

from rdkit import Chem
from collections import deque

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    
    The algorithm works by:
      1. Converting the SMILES to a molecule and adding explicit hydrogens.
      2. Identifying candidate carbon atoms (i.e. sp3, non‐ring carbons) that carry 
         at least one “free” or unmodified –OH group. Here we define an alcohol –OH 
         group as an oxygen atom (atomic number 8) that is attached to at least one 
         hydrogen and is not involved in a double bond (i.e. not a carbonyl oxygen).
      3. Building a graph among candidate carbons (an edge exists if two candidate 
         carbons are directly connected by a single, non‐ring bond).
      4. Scanning the candidate graph for a contiguous, acyclic path of at least 3 
         candidate carbons.
         
    Args:
      smiles (str): SMILES string of the molecule.
      
    Returns:
      (bool, str): A tuple where the first element is True if an alditol chain was found,
                   and False otherwise. The second element gives the reasoning.
    """
    # Convert SMILES to molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens – this will help in counting hydrogens on O atoms
    mol = Chem.AddHs(mol)
    
    # Define a helper function: is_alcohol.
    # We consider an oxygen as part of a free hydroxyl if:
    #   - It is O (atomic number 8)
    #   - It is not involved in any double bond (i.e. all its bonds are single)
    #   - It carries at least one hydrogen.
    def is_alcohol(oxygen_atom):
        if oxygen_atom.GetAtomicNum() != 8:
            return False
        # Check that none of its bonds is a double bond (avoids carbonyl)
        for bond in oxygen_atom.GetBonds():
            if bond.GetBondTypeAsDouble() == 2.0:
                return False
        # Check that it has at least one hydrogen (explicit hydrogens were added)
        if oxygen_atom.GetTotalNumHs() < 1:
            return False
        return True
    
    # Identify candidate carbons: sp3, non‐ring carbons that have at least one attached –OH group.
    candidate_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        # Only consider sp3 carbons that are not part of a ring.
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        if atom.IsInRing():
            continue
        # Check if any neighbor is an alcohol oxygen.
        has_alcohol = False
        for nb in atom.GetNeighbors():
            if is_alcohol(nb):
                has_alcohol = True
                break
        if has_alcohol:
            candidate_indices.append(atom.GetIdx())
    
    if not candidate_indices:
        return False, "No candidate carbon bearing a free –OH group found"
    
    # Build a graph connecting candidate carbons if they share a single bond.
    candidate_set = set(candidate_indices)
    graph = {idx: [] for idx in candidate_set}
    for idx in candidate_set:
        atom = mol.GetAtomWithIdx(idx)
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            # Connect only candidate carbons
            # Also, only consider bonds that are single and not in a ring.
            if nb_idx in candidate_set:
                bond = mol.GetBondBetweenAtoms(idx, nb_idx)
                if bond is not None and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                    # Occasionally, bonds may be flagged as in a ring;
                    # here we prefer to keep only acyclic connections.
                    if not bond.IsInRing():
                        if nb_idx not in graph[idx]:
                            graph[idx].append(nb_idx)
                        if idx not in graph[nb_idx]:
                            graph[nb_idx].append(idx)
    
    # A helper function to get connected components of an undirected graph.
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
    
    # Now we need to check each connected component for a linear (acyclic) path
    # of at least 3 candidate carbons.
    # We enumerate simple paths by DFS.
    def dfs_paths(g, current, visited):
        # yield the current path
        yield visited
        for nbr in g[current]:
            if nbr not in visited:
                for path in dfs_paths(g, nbr, visited + [nbr]):
                    yield path
    
    for comp in components:
        if len(comp) < 3:
            continue
        comp_graph = {node: graph[node] for node in comp}
        # For each candidate node in the component, try to find a path of length >=3.
        for start in comp_graph:
            for path in dfs_paths(comp_graph, start, [start]):
                if len(path) >= 3:
                    # All carbons along this path have already been preselected for having a free –OH.
                    return True, f"Found alditol chain of {len(path)} contiguous carbon(s) bearing free –OH groups."
    return False, "No linear alditol chain pattern detected"

# Example usage:
if __name__ == "__main__":
    # Test a few provided SMILES strings
    test_smiles = [
        "C([C@H]([C@H]([C@@H]([C@H]([13CH2]O)O)O)O)O)O",  # D-(1-(13)C)glucitol
        "OC[C@@H](O)[C@H](O)CO",                           # D-threitol
        "C[C@@](O)(CO)[C@@H](O)CO",                        # 2-methylerythritol
        "OC[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO",         # L-iditol
        "O([C@@H]([C@H](O[C@@H]1O[C@H]([C@@H](O)[C@@H](O)[C@@H]1O)C)[C@@H](O)CO)[C@H](O)CO)[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO",  # alpha-L-Fucp-(1->3)-[beta-D-Galp-(1->4)]-D-Glc-ol
        "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)C1OC(COC2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)CO)CO",  # 2-(D-glucopyranosyloxy)-3-hydroxypropyl L-fucopyranoside
        "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OCC(O)CO)CO",  # 1-O-(beta-D-glucosyl)glycerol
        "C(CO)(CO)O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO",  # 2-O-(beta-D-glucosyl)glycerol
        "OC[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO",         # D-iditol
        "OC[C@H](O)[C@H](O)[C@H](O)[C@@H](O)CO",          # L-altritol
        "OCC(O)CO"                                       # glycerol, smallest alditol
    ]
    for s in test_smiles:
        res, msg = is_alditol(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n{'-'*40}")