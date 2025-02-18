"""
Classifies: CHEBI:26244 prenols
"""
"""
Classifies: Prenols – any alcohol possessing the general formula 
    H-[CH2C(Me)=CHCH2]nOH
in which the carbon skeleton is exclusively composed of one or more isoprene units.

The following heuristic is used:
  1. The molecule must be valid, acyclic, and contain only C and O as heavy atoms.
  2. The prenol “backbone” is defined as the longest continuous chain (subgraph) of carbon atoms.
  3. The number of carbons in the backbone must be a multiple of 4 (each isoprene unit contributes 4 backbone C).
     (Thus, if backbone_length = L then the number isoprene units is n = L/4.)
  4. The number of carbon–carbon double bonds along the backbone must equal n.
  5. Exactly n methyl substituents (CH3 groups) must be attached to backbone carbons.
  6. The molecule must have exactly one free hydroxyl group ([OX2H]) that is directly attached to a terminal carbon of the backbone.
    
If all tests are passed, the molecule is classified as a prenol.
"""

from rdkit import Chem
from rdkit.Chem import rdchem
from collections import deque

def is_prenols(smiles: str):
    """
    Determines if a molecule is a prenol based on its SMILES.
    
    A prenol is defined as an alcohol whose carbon skeleton is exclusively
    composed of one or more isoprene units. Our algorithm:
      - Parses the molecule and verifies it is acyclic.
      - Requires heavy atoms to be only carbon and oxygen.
      - Extracts the longest chain from the carbon-only subgraph as the backbone.
      - Requires the backbone length be a multiple of 4 (n isoprene units, backbone length L = 4*n).
      - Requires that along that backbone, the number of C=C double bonds equals n.
      - Requires that exactly n CH3 substituents (methyl groups; defined as a carbon 
        having exactly 1 heavy neighbor) are attached off the backbone.
      - Requires exactly one free hydroxyl group ([OX2H]) and that its attached carbon is 
        one of the terminal carbons of the backbone.
    
    Args:
        smiles (str): SMILES string for the molecule.
    
    Returns:
        (bool, str): Tuple with True plus explanatory message if classified as prenol;
                     otherwise False with a reason.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # 1. The molecule must be acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule is cyclic, not a prenol"
        
    # 2. Molecule must contain only C and O as heavy atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (6, 8):
            return False, f"Molecule contains atom {atom.GetSymbol()} not in (C, O)"
     
    # 3. Create a subgraph of only carbon atoms.
    carbon_idxs = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_idxs:
        return False, "No carbon atoms found"
    
    # Build an adjacency dictionary for the carbon subgraph.
    # Key: carbon atom index; Value: list of neighboring carbon atom indices.
    carbon_adj = {idx: [] for idx in carbon_idxs}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            i1 = a1.GetIdx()
            i2 = a2.GetIdx()
            if i1 in carbon_adj and i2 in carbon_adj:
                carbon_adj[i1].append(i2)
                carbon_adj[i2].append(i1)
    
    # Identify terminal carbons (degree 1 in carbon subgraph).
    terminal_carbons = [idx for idx, nbrs in carbon_adj.items() if len(nbrs) == 1]
    if len(terminal_carbons) < 2:
        return False, "Could not identify two terminal carbons in the carbon skeleton"
    
    # 4. Find the longest path between any two terminals.
    # In acyclic graphs the unique simple path between any two nodes can be found by BFS.
    def bfs_path(start, goal):
        # Returns list of atom indices from start to goal.
        queue = deque([[start]])
        visited = set([start])
        while queue:
            path = queue.popleft()
            current = path[-1]
            if current == goal:
                return path
            for neighbor in carbon_adj[current]:
                if neighbor not in visited:
                    visited.add(neighbor)
                    new_path = path + [neighbor]
                    queue.append(new_path)
        return None

    longest_path = []
    for i in range(len(terminal_carbons)):
        for j in range(i+1, len(terminal_carbons)):
            start = terminal_carbons[i]
            end = terminal_carbons[j]
            path = bfs_path(start, end)
            if path and len(path) > len(longest_path):
                longest_path = path

    if not longest_path:
        return False, "Could not determine a longest carbon chain"
    
    # This longest path is our candidate backbone.
    backbone = longest_path
    backbone_length = len(backbone)
    
    # 5. For a prenol, the overall carbon count should be 5*n.
    # Since the true repeating unit is isoprene, the backbone (chain) consists of 4*n carbons;
    # the extra n are methyl groups attached to the backbone.
    if backbone_length % 4 != 0:
        return False, (f"Backbone length is {backbone_length}, which is not a multiple of 4; "
                       "cannot determine isoprene repeats")
    n_units = backbone_length // 4  # number of isoprene units expected
    
    # 6. Count C=C double bonds along the backbone.
    double_bonds_in_backbone = 0
    for i in range(len(backbone)-1):
        bond = mol.GetBondBetweenAtoms(backbone[i], backbone[i+1])
        if bond is not None and bond.GetBondType() == rdchem.BondType.DOUBLE:
            double_bonds_in_backbone += 1
    if double_bonds_in_backbone != n_units:
        return False, (f"Mismatch in double bond count along the backbone: found {double_bonds_in_backbone} "
                       f"double bond(s) but expected {n_units} (1 per isoprene repeat)")
    
    # 7. Count methyl substituents (CH3 groups) attached to the backbone.
    methyl_count = 0
    for idx in backbone:
        atom = mol.GetAtomWithIdx(idx)
        for nbr in atom.GetNeighbors():
            # Look only at carbon neighbors that are NOT part of the backbone.
            if nbr.GetAtomicNum() == 6 and nbr.GetIdx() not in backbone:
                # For a methyl group, the carbon should have only one heavy neighbor (the backbone carbon)
                heavy_neighbors = [a for a in nbr.GetNeighbors() if a.GetAtomicNum() > 1]
                if len(heavy_neighbors) == 1:
                    methyl_count += 1
    if methyl_count != n_units:
        return False, (f"Mismatch in methyl substituent count: found {methyl_count} "
                       f"but expected {n_units} (one per isoprene unit)")
    
    # 8. Identify free hydroxyl groups (using SMARTS [OX2H]).
    oh_pat = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pat)
    if not oh_matches:
        return False, "No free hydroxyl (–OH) group found"
    
    # In our prenol, we expect exactly one free –OH.
    if len(oh_matches) != 1:
        return False, (f"Expected exactly one free hydroxyl group, "
                       f"found {len(oh_matches)}")
    
    # 9. Ensure the –OH is attached to a terminal carbon of the backbone.
    oh_atom_idx = oh_matches[0][0]  # the oxygen atom index
    oh_atom = mol.GetAtomWithIdx(oh_atom_idx)
    # Its only heavy neighbor should be a carbon.
    oh_neighbors = [nbr for nbr in oh_atom.GetNeighbors() if nbr.GetAtomicNum() > 1]
    if not oh_neighbors:
        return False, "Hydroxyl group is isolated"
    attached_c = oh_neighbors[0]
    if attached_c.GetIdx() not in (backbone[0], backbone[-1]):
        return False, "Hydroxyl group is not attached to a terminal carbon of the backbone"
    
    # If all tests pass, report a successful classification.
    msg = (f"Classified as prenol: backbone has {backbone_length} carbons "
           f"({n_units} isoprene unit(s)), {double_bonds_in_backbone} C=C bond(s) along backbone, "
           f"{methyl_count} methyl substituent(s), and one terminal –OH group.")
    return True, msg

# For example usage (for testing):
if __name__ == "__main__":
    # Example: prenol ("prenol": CC(C)=CCO)
    test_smiles = "CC(C)=CCO"
    result, reason = is_prenols(test_smiles)
    print("Result:", result)
    print("Reason:", reason)