"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

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

    # Define the ester group pattern, R-C(=O)-O-R'
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    for match in ester_matches:
        # Get the index of the ester carbonyl carbon
        carbonyl_c_idx = match[0]

        # Use a BFS approach to find exactly 10 sequential carbons from carbonyl_c_idx
        def find_decanoic_chain(atom_idx):
            queue = [(atom_idx, 1)]  # Start with the first carbon (carbonyl carbon)
            visited = {atom_idx}
            while queue:
                current_idx, current_length = queue.pop(0)
                if current_length > 10:
                    continue
                if current_length == 10:
                    return True
                current_atom = mol.GetAtomWithIdx(current_idx)
                for neighbor in current_atom.GetNeighbors():
                    neighbor_idx = neighbor.GetIdx()
                    # Must be a carbon atom and not visited
                    if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited:
                        queue.append((neighbor_idx, current_length + 1))
                        visited.add(neighbor_idx)
            return False

        # Check both sides of the ester to find a 10-carbon acid chain
        for neighbor in mol.GetAtomWithIdx(carbonyl_c_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # carbon
                if find_decanoic_chain(neighbor.GetIdx()):
                    return True, "Contains ester group with a decanoic acid (10-carbon) chain"

    return False, "Ester group without a proper decanoic acid chain"