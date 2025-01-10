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

    # Iterate over all aliphatic nitrogen atoms (excluding aromatic and ring nitrogens)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic() and not atom.IsInRing():
            # Perform BFS to find a path to an aromatic ring via aliphatic carbons
            visited = set()
            queue = [(atom.GetIdx(), [])]  # Each element is (atom index, path)
            while queue:
                current_idx, path = queue.pop(0)
                if current_idx in visited:
                    continue
                visited.add(current_idx)
                current_atom = mol.GetAtomWithIdx(current_idx)

                if current_atom.GetIsAromatic() and current_atom.GetAtomicNum() == 6:
                    # Found an aromatic carbon
                    # Ensure the path consists only of aliphatic carbons
                    if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in path):
                        return True, "Contains amino group connected via alkyl chain to aromatic ring"
                else:
                    # Continue traversal through aliphatic carbons
                    if (current_atom.GetAtomicNum() == 6 and not current_atom.GetIsAromatic()) or current_idx == atom.GetIdx():
                        for neighbor in current_atom.GetNeighbors():
                            nbr_idx = neighbor.GetIdx()
                            if nbr_idx not in visited:
                                queue.append((nbr_idx, path + [current_idx]))
        else:
            continue

    # Check if nitrogen is directly attached to an aromatic ring (arylamine)
    arylamine_pattern = Chem.MolFromSmarts("[#7;!H0]-a")  # Nitrogen attached to aromatic atom
    if mol.HasSubstructMatch(arylamine_pattern):
        return False, "Amino group directly attached to aromatic ring, not an aralkylamine"

    return False, "Does not contain amino group connected via alkyl chain to aromatic ring"