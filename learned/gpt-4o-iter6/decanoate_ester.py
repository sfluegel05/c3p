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

    # Define the ester group pattern, e.g., R-C(=O)-O-R'
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Check for each ester if it can be classified as a decanoate ester
    for match in ester_matches:
        # Get the index of the ester carbonyl carbon
        carbonyl_c_idx = match[0]

        # Traverse from the ester carbonyl carbon to count carbons
        visited = {carbonyl_c_idx}
        def traverse(atom_idx, carbon_count):
            if carbon_count == 10:
                return True
            atom = mol.GetAtomWithIdx(atom_idx)
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                if neighbor_idx not in visited and neighbor.GetAtomicNum() == 6:  # carbon
                    visited.add(neighbor_idx)
                    if traverse(neighbor_idx, carbon_count + 1):
                        return True
                    visited.remove(neighbor_idx)
            return False

        for neighbor in mol.GetAtomWithIdx(carbonyl_c_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:
                visited.add(neighbor.GetIdx())
                if traverse(neighbor.GetIdx(), 1):
                    return True, "Contains ester group with a decanoic acid (10-carbon) chain"

    return False, "Ester group without a proper decanoic acid chain"