"""
Classifies: CHEBI:24026 fatty alcohol
"""
"""
Classifies: fatty alcohol
"""
from rdkit import Chem

def is_fatty_alcohol(smiles: str):
    """
    Determines if a molecule is a fatty alcohol based on its SMILES string.
    A fatty alcohol is an aliphatic alcohol consisting of a chain of 3 to greater than 27 carbon atoms.
    Fatty alcohols may be saturated or unsaturated and may be branched or unbranched.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty alcohol, False otherwise
        str: Reason for classification
    """
    from rdkit.Chem import rdmolops

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identify hydroxyl groups attached to sp3 carbons (exclude phenols)
    hydroxyl_pattern = Chem.MolFromSmarts("[C;H1,H2,H3;!$(C=[!#6])]-[OH]")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    if not hydroxyl_matches:
        return False, "No suitable hydroxyl groups attached to sp3 carbons"

    for match in hydroxyl_matches:
        carbon_idx = match[0]

        # Check if the carbon is in a ring
        carbon_atom = mol.GetAtomWithIdx(carbon_idx)
        # We now include ring carbons if they are part of an aliphatic chain

        # Use DFS to find the longest carbon chain starting from carbon_idx
        visited = set()
        max_chain_length = 0

        def dfs(atom_idx, length):
            nonlocal max_chain_length
            visited.add(atom_idx)
            atom = mol.GetAtomWithIdx(atom_idx)
            max_chain_length = max(max_chain_length, length)
            for neighbor in atom.GetNeighbors():
                nbr_idx = neighbor.GetIdx()
                nbr_atom = mol.GetAtomWithIdx(nbr_idx)
                if nbr_idx not in visited and nbr_atom.GetAtomicNum() == 6:
                    dfs(nbr_idx, length + 1)

        dfs(carbon_idx, 1)  # Start from length 1 (carbon attached to OH)

        if max_chain_length >= 3:
            # Additional check: Ensure the main chain is aliphatic (exclude aromatics)
            # Get the path of the longest chain
            path = []

            def dfs_path(atom_idx, current_path):
                current_path.append(atom_idx)
                atom = mol.GetAtomWithIdx(atom_idx)
                for neighbor in atom.GetNeighbors():
                    nbr_idx = neighbor.GetIdx()
                    nbr_atom = mol.GetAtomWithIdx(nbr_idx)
                    if nbr_atom.GetAtomicNum() == 6 and nbr_idx not in current_path:
                        dfs_path(nbr_idx, current_path.copy())

            dfs_path(carbon_idx, [])
            # Check if any atom in the path is aromatic
            aromatic = False
            for idx in path:
                if mol.GetAtomWithIdx(idx).GetIsAromatic():
                    aromatic = True
                    break
            if aromatic:
                continue  # Skip if any atom in the main chain is aromatic

            return True, f"Contains aliphatic alcohol with a chain of {max_chain_length} carbons"

    return False, "No aliphatic chain with sufficient length connected to hydroxyl group"