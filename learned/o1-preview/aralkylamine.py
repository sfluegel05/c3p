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

    # SMARTS pattern for aralkylamine:
    # Aliphatic nitrogen connected via any number (including zero) of non-aromatic carbons to an aromatic ring
    pattern = Chem.MolFromSmarts("[#7;!R;!H0][C;!R;!A]@[C;!R;!A]@[a]")

    # Try to match the pattern in the molecule
    if mol.HasSubstructMatch(pattern):
        return True, "Contains amino group connected via alkyl chain to aromatic ring"

    # Check for cases where nitrogen is directly connected to aromatic ring (should not be aralkylamine)
    direct_aryl_amino = Chem.MolFromSmarts("[#7]-[a]")
    if mol.HasSubstructMatch(direct_aryl_amino):
        return False, "Amino group directly attached to aromatic ring, not an aralkylamine"

    # As an additional check, find all aliphatic nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and not atom.GetIsAromatic():
            # Exclude nitrogens in rings
            if not atom.IsInRing():
                # Perform DFS to find aromatic ring through aliphatic carbons
                visited = set()
                stack = [(atom.GetIdx(), [])]
                while stack:
                    current_idx, path = stack.pop()
                    if current_idx in visited:
                        continue
                    visited.add(current_idx)
                    current_atom = mol.GetAtomWithIdx(current_idx)
                    if current_atom.GetIsAromatic():
                        # Reached an aromatic atom via path of aliphatic carbons
                        if all(mol.GetAtomWithIdx(idx).GetAtomicNum() == 6 and not mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in path):
                            return True, "Contains amino group connected via alkyl chain to aromatic ring"
                    elif current_atom.GetAtomicNum() == 6 and not current_atom.GetIsAromatic():
                        # Continue traversal through aliphatic carbons
                        neighbors = [nbr.GetIdx() for nbr in current_atom.GetNeighbors()]
                        for nbr_idx in neighbors:
                            if nbr_idx not in visited:
                                stack.append((nbr_idx, path + [current_idx]))
    # No aralkylamine substructure found
    return False, "Does not contain amino group connected via alkyl chain to aromatic ring"