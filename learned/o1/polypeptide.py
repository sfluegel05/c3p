"""
Classifies: CHEBI:15841 polypeptide
"""
"""
Classifies: polypeptide - A peptide containing ten or more amino acid residues.
"""

from rdkit import Chem
from collections import defaultdict, deque

def is_polypeptide(smiles: str):
    """
    Determines if a molecule is a polypeptide (a peptide containing ten or more amino acid residues)
    based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a polypeptide, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a more specific peptide bond SMARTS pattern
    # Peptide bond: N-C(=O)-C (amide nitrogen connected to carbonyl carbon connected to alpha carbon)
    peptide_bond_smarts = "[NX3][CX3](=O)[CX4]"

    peptide_bond_pattern = Chem.MolFromSmarts(peptide_bond_smarts)
    if peptide_bond_pattern is None:
        return False, "Invalid peptide bond SMARTS pattern"

    # Find peptide bonds in the molecule
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)

    # Build a graph of connected amino acid residues
    # Nodes are alpha carbon atoms (Cα)
    # Edges are peptide bonds between amino acids
    amino_acid_nodes = set()
    edges = []

    for match in peptide_bonds:
        n_idx = match[0]     # Nitrogen atom index
        c_o_idx = match[1]   # Carbonyl carbon atom index
        c_alpha_idx = match[2]  # Alpha carbon atom index

        # Add the alpha carbon as a node
        amino_acid_nodes.add(c_alpha_idx)

        # Get the nitrogen connected to this alpha carbon (may belong to next amino acid)
        c_alpha = mol.GetAtomWithIdx(c_alpha_idx)
        neighbor_nitrogens = [nbr for nbr in c_alpha.GetNeighbors() if nbr.GetAtomicNum() == 7 and nbr.GetIdx() != n_idx]
        for neighbor_n in neighbor_nitrogens:
            neighbor_n_idx = neighbor_n.GetIdx()
            # Check if the nitrogen belongs to another peptide bond
            if mol.HasSubstructMatch(Chem.MolFromSmarts(f"[#{neighbor_n_idx}][CX3](=O)[CX4]")):
                edges.append((c_alpha_idx, neighbor_n_idx))
                amino_acid_nodes.add(neighbor_n_idx)

    # Build the graph and find connected components
    graph = defaultdict(list)
    for a1, a2 in edges:
        graph[a1].append(a2)
        graph[a2].append(a1)

    # Find the largest connected component
    visited = set()
    max_component_size = 0
    for node in amino_acid_nodes:
        if node not in visited:
            # BFS to find all connected nodes
            queue = deque([node])
            visited.add(node)
            component_size = 1
            while queue:
                current = queue.popleft()
                for neighbor in graph[current]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
                        component_size += 1
            if component_size > max_component_size:
                max_component_size = component_size

    num_amino_acids = max_component_size

    if num_amino_acids >= 10:
        return True, f"Contains {num_amino_acids} amino acid residues connected via peptide bonds"
    else:
        return False, f"Contains {num_amino_acids} amino acid residues, less than 10 required for polypeptide"