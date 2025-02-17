"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: alditol

An alditol is an acyclic polyol having the general formula:
   HOCH2[CH(OH)]nCH2OH
(formally derivable from an aldose by reduction of its carbonyl group)

In the previous attempt the identification of candidate carbons depended on a helper function
(is_alcohol) that was too strict. Many of the SMILES examples were not finding any candidate carbon 
bearing a free –OH group because the oxygen test (checking for exactly one heavy neighbor and no double bonds)
failed due to explicit hydrogen handling or stereochemical annotations.
 
In this improved version we redefine is_alcohol to:
  • Check that the atom is oxygen (symbol "O"),
  • Verify that all bonds from the oxygen are single,
  • Count its neighbors: we require that it is attached to exactly one non‐hydrogen atom,
  • And that it is attached to at least one hydrogen (via an explicit neighbor).
  
Candidate carbons are those carbons (atomic number 6) that are sp³ hybridized and not in a ring,
and that have at least one neighboring oxygen that qualifies as a free hydroxyl.
Then we build a connectivity graph among these candidate carbons (linking directly bonded carbons via single bonds outside rings)
and search for any contiguous acyclic chain (via DFS) of at least three candidate carbons.
If such a chain exists, we classify the molecule as an alditol.
"""
from rdkit import Chem
from collections import deque

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    
    An alditol is defined as an acyclic polyol with the formula:
         HOCH2[CH(OH)]nCH2OH.
    We search for a chain of sp3 non‐ring carbons with one or more free hydroxyl groups.
    
    Args:
       smiles (str): SMILES string of the molecule.
       
    Returns:
       (bool, str): Tuple where the first element is True if an alditol chain is found,
                    and False otherwise. The second element provides the reasoning.
    """
    # Convert SMILES to molecule and add explicit hydrogens.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.AddHs(mol)
    
    # Helper function to decide if an oxygen atom is a free hydroxyl.
    # Conditions:
    #   • The atom must be oxygen.
    #   • All bonds from oxygen must be single.
    #   • It must be attached to exactly one non‐hydrogen atom.
    #   • It must be attached to at least one hydrogen.
    def is_alcohol(oxygen_atom):
        if oxygen_atom.GetSymbol() != "O":
            return False
        # Check all bonds are single bonds.
        for bond in oxygen_atom.GetBonds():
            if bond.GetBondType() != Chem.rdchem.BondType.SINGLE:
                return False
        # Count non-hydrogen neighbors.
        heavy_neighbors = [nb for nb in oxygen_atom.GetNeighbors() if nb.GetAtomicNum() != 1]
        if len(heavy_neighbors) != 1:
            return False
        # Count explicit hydrogen neighbors.
        hydrogen_neighbors = [nb for nb in oxygen_atom.GetNeighbors() if nb.GetAtomicNum() == 1]
        if len(hydrogen_neighbors) < 1:
            return False
        return True

    candidate_indices = []
    # Identify candidate carbon atoms:
    #   • Must be carbon (atomic number 6)
    #   • sp3 hybridized
    #   • Not in any ring (acyclic chain is desired)
    #   • Has at least one attached (free) hydroxyl group.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 6:
            continue
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        if atom.IsInRing():
            continue
        has_free_oh = False
        for nb in atom.GetNeighbors():
            # If neighbor qualifies as a free hydroxyl, mark candidate.
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
    
    # Get connected components from the candidate graph.
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
    
    # For each connected component, search for an acyclic chain of at least three candidate carbons.
    def dfs_paths(g, current, visited):
        yield visited
        for nbr in g[current]:
            if nbr not in visited:
                yield from dfs_paths(g, nbr, visited + [nbr])
    
    for comp in components:
        if len(comp) < 3:
            continue
        # Build subgraph for this component.
        comp_graph = {node: graph[node] for node in comp}
        for start in comp_graph:
            for path in dfs_paths(comp_graph, start, [start]):
                if len(path) >= 3:
                    return True, f"Found alditol chain of {len(path)} contiguous candidate carbons bearing free –OH groups."
    
    return False, "No linear alditol chain pattern detected"

# Example usage and testing: run this module to test several examples.
if __name__ == "__main__":
    test_smiles = [
        "C([C@H]([C@H]([C@@H]([C@H]([13CH2]O)O)O)O)O)O",  # D-(1-(13)C)glucitol
        "OC[C@@H](O)[C@H](O)CO",                           # D-threitol
        "C[C@@](O)(CO)[C@@H](O)CO",                        # 2-methylerythritol
        "OC[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)CO",         # L-iditol
        "O([C@@H]([C@H](O[C@@H]1O[C@H]([C@@H](O)[C@@H](O)[C@@H]1O)C)[C@@H](O)CO)[C@H](O)CO)[C@@H]2O[C@@H]([C@H](O)[C@H](O)[C@H]2O)CO",  # glycoside example
        "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)C1OC(COC2O[C@H]([C@@H](O)[C@@H](O)[C@@H]2O)C)CO)CO",  # another glycoside variant
        "O1[C@@H]([C@@H](O)[C@H](O)[C@@H](O)[C@@H]1OCC(O)CO)CO",  # 1-O-(beta-D-glucosyl)glycerol variant
        "C(CO)(CO)O[C@@H]1O[C@@H]([C@@H](O)[C@@H]([C@H]1O)O)CO",  # 2-O-(beta-D-glucosyl)glycerol variant
        "OC[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)CO",         # D-iditol
        "OC[C@H](O)[C@H](O)[C@H](O)[C@@H](O)CO",          # L-altritol
        "OCC(O)CO"                                       # glycerol
    ]
    for s in test_smiles:
        res, msg = is_alditol(s)
        print(f"SMILES: {s}\nResult: {res}\nReason: {msg}\n{'-'*60}")