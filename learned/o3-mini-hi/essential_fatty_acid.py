"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential fatty acid
Definition: A free (non-esterified) fatty acid that is acyclic, linear (unbranched), 
has exactly one terminal carboxylic acid group (acid carbon with exactly one carbon neighbor),
contains only C, H, and O (with exactly 2 O atoms, from the acid group),
has a sufficiently long chain (≥16 carbons), and has at least 2 carbon–carbon double bonds aside from the acid carbonyl.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is a free essential fatty acid based on its SMILES string.
    Criteria:
      - Valid SMILES.
      - Acyclic (no rings).
      - Contains only C, H, and O (extra heteroatoms not allowed).
      - Contains exactly one free carboxylic acid group. This group must be terminal:
        the carboxyl carbon should have exactly one carbon neighbor.
      - The molecule must be linear (i.e. unbranched): the longest chain of carbons should
        include all carbons present.
      - Contains at least 16 carbon atoms.
      - Excluding the acid carbonyl, the molecule must have at least 2 carbon–carbon double bonds.
      
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule meets the criteria, False otherwise.
        str: Explanation for the decision.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check that the molecule is acyclic.
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings, not a simple acyclic fatty acid"
    
    # Ensure only allowed elements: hydrogen (1), carbon (6), and oxygen (8)
    allowed_atomic_nums = {1, 6, 8}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in allowed_atomic_nums:
            return False, "Molecule has additional heteroatoms (e.g., N, P) not found in a simple fatty acid"
    
    # Look for the free carboxylic acid group.
    # SMARTS matches [CX3](=O)[O;H1,-] so it captures both protonated and ionized forms.
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    acid_carbon_indices = set(match[0] for match in acid_matches)
    if len(acid_carbon_indices) != 1:
        return False, f"Expected one free carboxylic acid group, found {len(acid_carbon_indices)}"
    
    acid_carbon_idx = next(iter(acid_carbon_indices))
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    # Check that the acid carbon is terminal (has exactly one carbon neighbor)
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxylic acid group is not terminal (acid carbon should have exactly one carbon neighbor)"
    
    # Create a carbon-only graph to assess linearity.
    # Build a dictionary of carbon indices mapping to neighbors (only carbons).
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_indices:
        return False, "No carbon atoms found"
    carbon_graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            # Add neighbors if both atoms are in our carbon graph.
            if idx1 in carbon_graph and idx2 in carbon_graph:
                carbon_graph[idx1].append(idx2)
                carbon_graph[idx2].append(idx1)
    
    # For an acyclic graph (tree) the longest path can be found by:
    # (1) pick an arbitrary carbon and perform BFS to find the farthest carbon.
    # (2) from that farthest, perform BFS again to get the longest distance.
    def bfs_farthest(start, graph):
        visited = {start}
        queue = [(start, 0)]
        farthest_node = start
        max_dist = 0
        while queue:
            current, dist = queue.pop(0)
            if dist > max_dist:
                max_dist = dist
                farthest_node = current
            for nbr in graph[current]:
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, dist+1))
        return farthest_node, max_dist
    
    arbitrary = carbon_indices[0]
    node1, _ = bfs_farthest(arbitrary, carbon_graph)
    node2, longest_chain_len = bfs_farthest(node1, carbon_graph)
    
    total_carbons = len(carbon_indices)
    # If the longest chain (number of bonds+1) is not equal to the total number of carbons,
    # then there is branching.
    if longest_chain_len + 1 != total_carbons:
        return False, "Molecule is branched; not a simple linear fatty acid"
    
    # Check total number of carbons meets chain length criteria.
    if total_carbons < 16:
        return False, f"Not enough carbon atoms (found {total_carbons}, need at least 16)"
    
    # Count C–C double bonds that are not part of the acid carbonyl
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # exclude double bonds involving the acid carbon
                if a1.GetIdx() == acid_carbon_idx or a2.GetIdx() == acid_carbon_idx:
                    continue
                cc_double_bonds += 1
    if cc_double_bonds < 2:
        return False, f"Not polyunsaturated enough (found {cc_double_bonds} C–C double bond(s); need at least 2)"
    
    return True, (f"Contains a terminal free carboxylic acid group, {total_carbons} carbons, and {cc_double_bonds} "
                  "carbon–carbon double bonds, consistent with an essential fatty acid")

# Example usage:
if __name__ == "__main__":
    test_smiles = "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC(O)=O"  # all-cis-octadeca-6,9,12,15-tetraenoic acid
    result, reason = is_essential_fatty_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)