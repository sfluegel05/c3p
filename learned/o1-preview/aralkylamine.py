"""
Classifies: CHEBI:18000 aralkylamine
"""
"""
Classifies: aralkylamine
"""

from rdkit import Chem

def is_aralkylamine(smiles: str):
    """
    Determines if a molecule is an aralkylamine based on its SMILES string.
    An aralkylamine is defined as an alkylamine in which the alkyl group is substituted by an aromatic group.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an aralkylamine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Iterate over all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Check if nitrogen is aliphatic (not aromatic)
            if not atom.GetIsAromatic():
                # Check if nitrogen is not in a ring
                if not atom.IsInRing():
                    # Perform BFS traversal to find aromatic ring via aliphatic carbons
                    visited = set()
                    queue = []
                    atom_idx = atom.GetIdx()
                    queue.append((atom_idx, []))  # (current atom idx, path)

                    while queue:
                        current_idx, path = queue.pop(0)
                        if current_idx in visited:
                            continue
                        visited.add(current_idx)
                        current_atom = mol.GetAtomWithIdx(current_idx)
                        # Skip the nitrogen atom itself in path
                        if current_idx != atom_idx:
                            # Check if current atom is aromatic
                            if current_atom.GetIsAromatic():
                                # Check that path consists only of aliphatic carbons
                                is_valid_path = True
                                for idx in path:
                                    path_atom = mol.GetAtomWithIdx(idx)
                                    if path_atom.GetAtomicNum() != 6:
                                        is_valid_path = False
                                        break
                                    if path_atom.GetIsAromatic():
                                        is_valid_path = False
                                        break
                                if is_valid_path and len(path) >= 1:
                                    return True, "Contains amino group connected via alkyl chain to aromatic ring"
                            # Continue traversal only through aliphatic carbons
                            elif current_atom.GetAtomicNum() == 6 and not current_atom.GetIsAromatic():
                                # Visit neighbors
                                for neighbor in current_atom.GetNeighbors():
                                    neighbor_idx = neighbor.GetIdx()
                                    if neighbor_idx not in visited:
                                        queue.append((neighbor_idx, path + [current_idx]))
                    # No valid path found for this nitrogen atom
    # No aralkylamine substructure found
    return False, "Does not contain amino group connected via alkyl chain to aromatic ring"