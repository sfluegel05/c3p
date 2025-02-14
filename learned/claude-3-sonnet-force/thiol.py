"""
Classifies: CHEBI:29256 thiol
"""
"""
Classifies: CHEBI:34831 thiol
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_thiol(smiles: str) -> tuple[bool, str]:
    """
    Determines if a molecule is a thiol based on its SMILES string.
    A thiol is an organosulfur compound in which a thiol group (-SH) is attached to a carbon atom
    of any aliphatic or aromatic moiety.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a thiol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for thiol group (-SH) pattern
    thiol_pattern = Chem.MolFromSmarts("[SH]")
    thiol_matches = mol.GetSubstructMatches(thiol_pattern)

    # Check if there is at least one thiol group
    if len(thiol_matches) == 0:
        return False, "No thiol group (-SH) found"

    def traverse_ring_or_chain(atom, visited):
        """
        Recursive function to traverse a ring or aliphatic chain starting from a given atom.
        Returns True if the traversal encounters a non-carbon atom or reaches a termination condition.
        """
        if atom.GetAtomicNum() != 6:  # Not a carbon atom
            return True

        visited.add(atom.GetIdx())
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited:
                if traverse_ring_or_chain(neighbor, visited):
                    return True

        return False

    for thiol_idx in thiol_matches:
        thiol_atom = mol.GetAtomWithIdx(thiol_idx)
        for neighbor_atom in thiol_atom.GetNeighbors():
            if neighbor_atom.GetAtomicNum() == 6:  # Carbon
                visited = set()
                if traverse_ring_or_chain(neighbor_atom, visited):
                    return True, "Contains a thiol group (-SH) attached to a carbon atom of an aliphatic or aromatic moiety"

    return False, "Thiol group not attached to a carbon atom of an aliphatic or aromatic moiety"