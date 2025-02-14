"""
Classifies: CHEBI:10036 wax ester
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_wax_ester(smiles: str):
    """
    Determines if a molecule is a wax ester based on its SMILES string.
    A wax ester is a fatty acid ester resulting from the condensation of the carboxy group of
    a fatty acid with the alcoholic hydroxy group of a fatty alcohol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a wax ester, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for ester linkage (-C(=O)O-)
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"

    # Function to count contiguous linear carbon atoms.
    def count_linear_carbon_chain(atom_idx, visited=set()):
        stack = [atom_idx]
        carbon_count = 0
        is_linear = True

        while stack and is_linear:
            atom_idx = stack.pop()
            if atom_idx in visited:
                is_linear = False
                continue
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:  # Continue if not a carbon atom
                continue
            carbon_count += 1
            neighbors = atom.GetNeighbors()
            
            if len(neighbors) > 2:
                is_linear = False  # Not linear if more than 2 neighbors

            for neighbor in neighbors:
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    stack.append(neighbor_idx)

        return carbon_count, is_linear

    # Verify that each ester group connects two sufficiently long and linear carbon chains
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    for match in ester_matches:
        atom_idx_c1 = match[0]  # Carbon of the carbonyl group
        atom_idx_o = match[2]   # Oxygen in -C(=O)O-

        # Get neighbors excluding the ones directly involved in the ester linkage
        c1_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(atom_idx_c1).GetNeighbors() if n.GetIdx() != match[1]]
        o_neighbors = [n.GetIdx() for n in mol.GetAtomWithIdx(atom_idx_o).GetNeighbors() if n.GetIdx() != match[1]]

        # Check if neighbors exist
        if not c1_neighbors or not o_neighbors:
            return False, "Missing chains connected to the ester linkage"

        c1_chain_lengths = [count_linear_carbon_chain(neigh_idx) for neigh_idx in c1_neighbors]
        o_chain_lengths = [count_linear_carbon_chain(neigh_idx) for neigh_idx in o_neighbors]

        # Find the longest linear chains connected to the ester linkage
        longest_c1_chain = max(c1_chain_lengths, key=lambda x: x[0])
        longest_o_chain = max(o_chain_lengths, key=lambda x: x[0])

        # Ensure both chains are sufficiently long (e.g., >=8 carbons) and linear
        if longest_c1_chain[0] < 8 or longest_o_chain[0] < 8 or not longest_c1_chain[1] or not longest_o_chain[1]:
            return False, f"Chains attached to ester group are too short or not linear (lengths found: {longest_c1_chain[0]} and {longest_o_chain[0]})"

    return True, "Contains an ester linkage with sufficiently long and linear carbon chains typical of wax esters"