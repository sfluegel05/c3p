"""
Classifies: CHEBI:17135 long-chain fatty alcohol
"""
from rdkit import Chem
from rdkit.Chem import rdmolops

def is_long_chain_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a long-chain fatty alcohol based on its SMILES string.
    A long-chain fatty alcohol is a fatty alcohol with a chain length ranging from C13 to C22.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a long-chain fatty alcohol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for presence of hydroxyl group (-OH)
    oh_group = Chem.MolFromSmarts('[OX2H]')
    if not mol.HasSubstructMatch(oh_group):
        return False, "No hydroxyl group (-OH) found"

    # Build a list of carbon atoms
    carbon_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]

    # Build a graph of carbon atoms connected via single bonds (excluding ring bonds)
    carbon_bonds = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        if bond.IsInRing():
            continue  # Exclude ring bonds
        if a1.GetAtomicNum() == 6 and a2.GetAtomicNum() == 6:  # Both atoms are carbon
            carbon_bonds.append((a1.GetIdx(), a2.GetIdx()))

    # Create a graph of carbon atoms
    import networkx as nx
    G = nx.Graph()
    G.add_nodes_from(carbon_atoms)
    G.add_edges_from(carbon_bonds)

    # Find the longest path in the carbon atom graph
    # For trees (acyclic graphs), we can find it using BFS
    def longest_path_length(graph):
        max_length = 0
        for node in graph.nodes():
            lengths = nx.single_source_shortest_path_length(graph, node)
            if lengths:
                current_max = max(lengths.values())
                if current_max > max_length:
                    max_length = current_max
        return max_length + 1  # Number of nodes in the longest path

    if G.number_of_nodes() == 0:
        return False, "No carbon chain found"

    longest_chain_length = longest_path_length(G)

    if longest_chain_length < 13:
        return False, f"Longest carbon chain length is {longest_chain_length}, less than 13"
    if longest_chain_length > 22:
        return False, f"Longest carbon chain length is {longest_chain_length}, greater than 22"

    # Check if hydroxyl group is attached to the chain
    # Get indices of atoms in the longest chain
    paths = dict(nx.all_pairs_shortest_path(G))
    longest_path_nodes = set()
    max_path_length = 0
    for source, target_dict in paths.items():
        for target, path in target_dict.items():
            if len(path) > max_path_length:
                max_path_length = len(path)
                longest_path_nodes = set(path)
    # Get the hydroxyl group attachment points
    oh_matches = mol.GetSubstructMatches(oh_group)
    attached = False
    for match in oh_matches:
        oxygen_idx = match[0]
        oxygen_atom = mol.GetAtomWithIdx(oxygen_idx)
        for neighbor in oxygen_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() in longest_path_nodes:
                attached = True
                break
        if attached:
            break
    if not attached:
        return False, "Hydroxyl group is not attached to the longest carbon chain"

    return True, "Molecule is a long-chain fatty alcohol with chain length between 13 and 22 carbons"