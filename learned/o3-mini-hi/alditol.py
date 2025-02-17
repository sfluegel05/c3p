"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol

An alditol is an acyclic polyol having the general formula:
  HOCH2[CH(OH)]nCH2OH
(i.e. it is formally derivable from an aldose by reduction of its carbonyl group).

This implementation looks for a contiguous, acyclic carbon chain (of at least three carbon atoms)
in which every carbon bears at least one “free” hydroxyl group. A free hydroxyl (–OH) is defined as:
  • An oxygen atom (atomic number 8) that, after adding explicit hydrogens, has a degree of 2 
    (i.e. the only heavy‐atom neighbor is the carbon to which it is bound)
  • And that has at least one hydrogen (and no extra heavy–atom neighbors).

Then among that candidate chain we require that the terminal (end) carbons behave “primarily”
(ideally CH2–OH; we require at least one hydrogen on the candidate carbon, but preferably two)
and that the interior carbons are “secondary” (ideally CH(OH); here we expect the number of hydrogens
to be no more than one). (Note that many “glycosidated” molecules will fail these tests, and the routine
may return (False, ...) for them.) 

The idea is to avoid selecting fragments from esters, phosphates or other modified –OH groups.
"""

from rdkit import Chem
from collections import deque

def get_connected_components(graph):
    """
    Given a graph as a dictionary {node: [neighbors]}, returns a list of sets,
    each representing a connected component.
    """
    seen = set()
    comps = []
    for node in graph:
        if node not in seen:
            comp = set()
            queue = deque([node])
            seen.add(node)
            while queue:
                current = queue.popleft()
                comp.add(current)
                for nbr in graph[current]:
                    if nbr not in seen:
                        seen.add(nbr)
                        queue.append(nbr)
            comps.append(comp)
    return comps

def is_free_oh(oxygen):
    """
    Determines if an oxygen atom is a free hydroxyl group.
    Criteria:
      - oxygen has atomic num = 8
      - its degree (number of neighbors) is exactly 2
      - exactly one neighbor is a hydrogen (atomic num 1)
      - (thus the other neighbor is the candidate carbon)
    """
    if oxygen.GetAtomicNum() != 8:
        return False
    # Use explicit hydrogens (we already added them)
    if oxygen.GetDegree() != 2:
        return False
    h_neighbors = sum(1 for nb in oxygen.GetNeighbors() if nb.GetAtomicNum() == 1)
    # Should have exactly one hydrogen neighbor
    if h_neighbors != 1:
        return False
    # Also make sure that among the heavy atoms it is attached only to one atom.
    heavy_neighbors = [nb for nb in oxygen.GetNeighbors() if nb.GetAtomicNum() != 1]
    if len(heavy_neighbors) != 1:
        return False
    return True

def enumerate_paths(graph, start, path=None, visited=None, max_length=10):
    """
    Recursively enumerate all simple paths (without repeated nodes)
    in the candidate graph starting from 'start' up to max_length.
    """
    if path is None:
        path = [start]
    if visited is None:
        visited = set([start])
    yield list(path)
    if len(path) >= max_length:
        return
    for nbr in graph[start]:
        if nbr in visited:
            continue
        # extend path
        new_path = list(path)
        new_path.append(nbr)
        new_visited = set(visited)
        new_visited.add(nbr)
        for p in enumerate_paths(graph, nbr, new_path, new_visited, max_length):
            yield p

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    
    We add explicit hydrogens and then search for a contiguous, acyclic chain of carbon
    atoms (of length at least 3) in which every carbon bears at least one free (non‐esterified, 
    non‐phosphorylated) hydroxyl group.
    
    Terminal carbons are required to have at least one hydrogen (ideally 2 for CH2–OH);
    interior carbons are expected to have at most one hydrogen (CHOH).
    
    Args:
      smiles (str): the SMILES string of the molecule.
      
    Returns:
      (bool, str): Classification result and reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    # Add explicit hydrogens
    mol = Chem.AddHs(mol)

    # First, mark candidate carbons: those bearing at least one free OH group.
    candidate_indices = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:  # consider only carbons
            continue
        # Check neighbors for free hydroxyl (–OH)
        has_free_oh = False
        for nb in atom.GetNeighbors():
            if nb.GetAtomicNum() == 8 and is_free_oh(nb):
                has_free_oh = True
                break
        if has_free_oh:
            candidate_indices.append(atom.GetIdx())
    if not candidate_indices:
        return False, "No carbon with a free –OH group detected"

    # Build a graph among candidate carbon atoms (if they share a bond).
    candidate_set = set(candidate_indices)
    graph = {idx: [] for idx in candidate_set}
    for idx in candidate_set:
        atom = mol.GetAtomWithIdx(idx)
        for nb in atom.GetNeighbors():
            nb_idx = nb.GetIdx()
            if nb_idx in candidate_set:
                graph[idx].append(nb_idx)

    # Get connected components in the candidate graph.
    components = get_connected_components(graph)
    # We require the alditol chain to have at least 3 carbons (as in HOCH2-CHOH-CH2OH)
    for comp in components:
        if len(comp) < 3:
            continue
        # Within this component we search for a linear (acyclic) path that might display
        # the canonical alditol pattern. Because our candidate graph (by definition) may be branched,
        # we enumerate all simple paths within the component (limiting path length to say 10).
        # We then check the free –OH and hydrogen counts on each carbon.
        comp_list = list(comp)
        comp_graph = {node: graph[node] for node in comp}
        found_path = None
        for start in comp_list:
            for path in enumerate_paths(comp_graph, start, max_length= len(comp)):
                if len(path) < 3:
                    continue
                # Make sure the atoms in the path are not part of a ring (the chain should be open).
                if any(mol.GetAtomWithIdx(idx).IsInRing() for idx in path):
                    continue
                valid_chain = True
                # Check each carbon in the chain: it must have at least one free –OH.
                # Also, enforce hydrogen count heuristics:
                #   For terminal atoms: require at least 1 hydrogen (ideally 2 are expected for CH2–OH)
                #   For interior atoms: require at most 1 hydrogen (as in CH(OH)) 
                for i, idx in enumerate(path):
                    atom = mol.GetAtomWithIdx(idx)
                    # check free –OH exists (we already required that for candidacy, but recheck if needed)
                    has_oh = any(nb.GetAtomicNum() == 8 and is_free_oh(nb) for nb in atom.GetNeighbors())
                    if not has_oh:
                        valid_chain = False
                        break
                    h_count = atom.GetTotalNumHs()  # after adding explicit hydrogens
                    if i == 0 or i == len(path)-1:
                        # Terminal carbon – expect primary alcohol behavior (ideally CH2–OH)
                        if h_count < 1:  # require at least one hydrogen (most cases show 2)
                            valid_chain = False
                            break
                    else:
                        # Interior carbon – expect secondary alcohol (ideally CH(OH))
                        if h_count > 1:
                            valid_chain = False
                            break
                if valid_chain:
                    found_path = path
                    break
            if found_path:
                break
        if found_path:
            n = len(found_path)
            return True, f"Found alditol chain with {n} contiguous carbon atoms bearing free –OH groups (terminal and interior carbons match expected hydrogen counts)."
    return False, "No linear alditol chain pattern detected"

# Example usage:
if __name__ == "__main__":
    # Some examples – you can add more SMILES strings to test the algorithm.
    test_smiles_list = [
        "C([C@H]([C@H]([C@@H]([C@H]([13CH2]O)O)O)O)O)O",  # D-(1-(13)C)glucitol (TP)
        "OC[C@@H](O)[C@H](O)CO",                          # D-threitol (TP)
        "C[C@@](O)(CO)[C@@H](O)CO",                       # 2-methylerythritol (was FN before)
        "OC(=O)CCCCCCCCCCCCCCC.OC(CO)CO",                 # Glycerol 1-hexadecanoate (FP, should be rejected)
        "OCCO"                                          # ethylene glycol (FP, too short)
    ]
    for s in test_smiles_list:
        res, msg = is_alditol(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n{'-'*40}")