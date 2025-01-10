"""
Classifies: CHEBI:17984 acyl-CoA
"""
"""
Classifies: acyl-CoA
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_acyl_CoA(smiles: str):
    """
    Determines if a molecule is an acyl-CoA based on its SMILES string.
    An acyl-CoA is a thioester resulting from the condensation of the thiol group of coenzyme A
    with the carboxy group of any carboxylic acid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an acyl-CoA, False otherwise
        str: Reason for classification
    """

    # Parse SMILES of input molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define thioester SMARTS pattern (C(=O)S)
    thioester_smarts = '[CX3](=O)[SX2]'
    thioester_pattern = Chem.MolFromSmarts(thioester_smarts)

    # Search for thioester groups
    thioester_matches = mol.GetSubstructMatches(thioester_pattern, useChirality=False)
    if not thioester_matches:
        return False, "No thioester group found"

    # Define adenine ring SMARTS pattern
    adenine_smarts = 'c1ncnc2ncnc12'
    adenine_pattern = Chem.MolFromSmarts(adenine_smarts)
    if adenine_pattern is None:
        return False, "Error in creating adenine pattern"

    # Search for adenine ring
    adenine_matches = mol.GetSubstructMatches(adenine_pattern, useChirality=False)
    if not adenine_matches:
        return False, "Adenine moiety not found"

    # Get atom indices of adenine ring
    adenine_atom_indices = set()
    for match in adenine_matches:
        adenine_atom_indices.update(match)

    # For each thioester group, check connectivity to adenine
    for thioester_match in thioester_matches:
        sulfur_idx = thioester_match[2]  # Index of sulfur atom in thioester group
        sulfur_atom = mol.GetAtomWithIdx(sulfur_idx)

        # Perform BFS from sulfur atom to check for path to adenine atoms
        visited = set()
        queue = [sulfur_atom.GetIdx()]
        found_path_to_adenine = False

        while queue:
            current_idx = queue.pop(0)
            if current_idx in visited:
                continue
            visited.add(current_idx)

            if current_idx in adenine_atom_indices:
                found_path_to_adenine = True
                break

            current_atom = mol.GetAtomWithIdx(current_idx)
            neighbors = current_atom.GetNeighbors()
            for neighbor in neighbors:
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited:
                    queue.append(neighbor_idx)

        if found_path_to_adenine:
            return True, "Contains thioester linked to CoA moiety (adenine ring connected via sulfur)"

    return False, "Thioester group not connected to CoA moiety"