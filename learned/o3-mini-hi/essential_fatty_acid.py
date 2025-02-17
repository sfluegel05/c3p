"""
Classifies: CHEBI:59549 essential fatty acid
"""
"""
Classifies: Essential fatty acid
Definition: A free (non-esterified) fatty acid that is acyclic, linear (unbranched),
contains only C, H, and O (with exactly 2 O atoms coming only from the terminal carboxyl group),
has exactly one terminal free carboxylic acid group (the acid carbon has exactly one carbon neighbor),
has a sufficiently long linear chain (>=16 carbons),
and (aside from the acid carbonyl) contains at least 2 carbon–carbon double bonds.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_essential_fatty_acid(smiles: str):
    """
    Determines if a molecule is a free essential fatty acid based on its SMILES string.
    Criteria:
      - Valid SMILES.
      - Acyclic (no rings).
      - Contains only C, H, and O (and exactly 2 O atoms overall)
      - Contains exactly one terminal free carboxylic acid group;
        the acid carbon must have exactly one carbon neighbor.
      - Is linear (i.e. unbranched): the longest carbon chain must account for all carbons.
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
    
    # Check allowed elements: only C (6), H (1), and O (8)
    allowed_atomic_nums = {1, 6, 8}
    oxygens = 0
    for atom in mol.GetAtoms():
        an = atom.GetAtomicNum()
        if an not in allowed_atomic_nums:
            return False, "Molecule has additional heteroatoms (e.g., N, P) not found in a simple fatty acid"
        if an == 8:
            oxygens += 1

    # For a free fatty acid, there should be exactly 2 oxygen atoms (from the COOH group)
    if oxygens != 2:
        return False, f"Molecule has {oxygens} oxygen atoms; a free fatty acid should have exactly 2 (from the carboxyl group)"

    # Find the free carboxylic acid group.
    # This SMARTS pattern captures a carboxyl group showing C(=O)[O;H1,-]
    acid_smarts = Chem.MolFromSmarts("[CX3](=O)[O;H1,-]")
    acid_matches = mol.GetSubstructMatches(acid_smarts)
    acid_carbon_indices = set(match[0] for match in acid_matches)
    if len(acid_carbon_indices) != 1:
        return False, f"Expected one free carboxylic acid group, found {len(acid_carbon_indices)}"
    
    acid_carbon_idx = next(iter(acid_carbon_indices))
    acid_carbon = mol.GetAtomWithIdx(acid_carbon_idx)
    # Verify that the acid carbon is terminal (should have exactly one carbon neighbor)
    carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
    if len(carbon_neighbors) != 1:
        return False, "The carboxylic acid group is not terminal (acid carbon should have exactly one carbon neighbor)"
    
    # Assess linearity: Build a graph of only carbon atoms
    carbon_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if not carbon_indices:
        return False, "No carbon atoms found in the molecule"
    
    # Build a dictionary representing the carbon graph: each node connected to carbon neighbors only.
    carbon_graph = {idx: [] for idx in carbon_indices}
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
            idx1 = a1.GetIdx()
            idx2 = a2.GetIdx()
            if idx1 in carbon_graph and idx2 in carbon_graph:
                carbon_graph[idx1].append(idx2)
                carbon_graph[idx2].append(idx1)
    
    # Use breadth-first search (BFS) to find the longest path in the tree of carbons.
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
                    queue.append((nbr, dist + 1))
        return farthest_node, max_dist

    arbitrary = carbon_indices[0]
    node1, _ = bfs_farthest(arbitrary, carbon_graph)
    node2, longest_path_length = bfs_farthest(node1, carbon_graph)
    total_carbons = len(carbon_indices)
    # For a linear chain, the longest path (in terms of number of bonds) plus one should equal the total number of carbons.
    if longest_path_length + 1 != total_carbons:
        return False, "Molecule is branched; not a simple linear fatty acid"

    # Check that the total number of carbon atoms meets the length criteria.
    if total_carbons < 16:
        return False, f"Not enough carbon atoms (found {total_carbons}, need at least 16)"

    # Count carbon-carbon double bonds that are not part of the acid carbonyl.
    cc_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            a1 = bond.GetBeginAtom()
            a2 = bond.GetEndAtom()
            if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:
                # Exclude bonds that involve the acid carbon.
                if a1.GetIdx() == acid_carbon_idx or a2.GetIdx() == acid_carbon_idx:
                    continue
                cc_double_bonds += 1
    if cc_double_bonds < 2:
        return False, f"Not polyunsaturated enough (found {cc_double_bonds} C–C double bond(s); need at least 2)"

    return True, (f"Contains a terminal free carboxylic acid group, {total_carbons} carbons, and {cc_double_bonds} "
                  "carbon–carbon double bonds, consistent with an essential fatty acid")

# Example usage:
if __name__ == "__main__":
    # Test with one of the provided examples:
    test_smiles = "CC\\C=C/C\\C=C/C\\C=C/C\\C=C/CCCCC(O)=O"  # all-cis-octadeca-6,9,12,15-tetraenoic acid
    result, reason = is_essential_fatty_acid(test_smiles)
    print("Result:", result)
    print("Reason:", reason)