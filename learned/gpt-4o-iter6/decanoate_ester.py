"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem
from rdkit.Chem import rdqueries

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has an ester functional group with a decanoic acid chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester group pattern
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Check for each ester if it can be classified as a decanoate ester
    for match in ester_matches:
        # Get the indices of the ester functional group
        carbonyl_c, ester_o = match[0], match[2]

        # Explore chains emerging from oxygen to find decanoic carbon chain
        atom = mol.GetAtomWithIdx(ester_o)
        visited_atoms = {ester_o}
        
        def traverse(atom_idx, depth):
            if depth == 10:
                return True
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() != 6:  # Check if current atom is carbon
                return False
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited_atoms:
                    visited_atoms.add(neighbor_idx)
                    if traverse(neighbor_idx, depth+1):
                        return True
                    visited_atoms.remove(neighbor_idx)
            return False
        
        # Start traversal from the atom next to ester oxygen
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Must be Carbon
                visited_atoms.add(neighbor.GetIdx())
                if traverse(neighbor.GetIdx(), 1):
                    return True, "Contains ester group with a decanoic acid (10-carbon) chain"

    return False, "Ester group without a proper decanoic acid chain"