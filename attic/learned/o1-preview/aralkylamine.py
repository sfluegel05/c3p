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

    # Iterate over all nitrogen atoms (excluding aromatic nitrogens)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic():
            # Skip if nitrogen is directly connected to an aromatic carbon (arylamine)
            if any(neighbor.GetIsAromatic() and neighbor.GetAtomicNum() == 6 for neighbor in atom.GetNeighbors()):
                continue

            # Perform BFS to find a path to an aromatic ring via non-aromatic atoms
            visited = set()
            queue = [(atom.GetIdx(), [atom.GetIdx()])]  # Each element is (atom index, path)
            while queue:
                current_idx, path = queue.pop(0)
                if current_idx in visited:
                    continue
                visited.add(current_idx)
                current_atom = mol.GetAtomWithIdx(current_idx)

                if current_idx != atom.GetIdx() and current_atom.GetIsAromatic() and current_atom.GetAtomicNum() == 6:
                    # Found an aromatic carbon
                    return True, "Contains amino group connected via alkyl chain to aromatic ring"
                else:
                    # Continue traversal through non-aromatic atoms
                    if not current_atom.GetIsAromatic():
                        for neighbor in current_atom.GetNeighbors():
                            nbr_idx = neighbor.GetIdx()
                            if nbr_idx not in visited:
                                queue.append((nbr_idx, path + [current_idx]))
    return False, "Does not contain amino group connected via alkyl chain to aromatic ring"