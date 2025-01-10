"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester has an ester functional group with a decanoic acid side chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES into an RDKit molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the ester functional group pattern: R-C(=O)-O-R'
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # Function to find a consecutive chain of specified carbons starting from a given atom
    def find_carbon_chain(atom_idx, length=10):
        visited = set()
        stack = [(atom_idx, 0)]

        while stack:
            current_idx, count = stack.pop()
            if count == length:
                return True
            visited.add(current_idx)
            for neighbor in mol.GetAtomWithIdx(current_idx).GetNeighbors():
                neighbor_idx = neighbor.GetIdx()
                # Ensure we continue on carbon atoms and not visited previously
                if neighbor.GetAtomicNum() == 6 and neighbor_idx not in visited:
                    stack.append((neighbor_idx, count + 1))
        return False

    # Check ester group occurrences for decanoic chain
    for match in ester_matches:
        carbonyl_c_idx = match[0]  # The ester carbonyl carbon
        oxygen_idx = match[2]  # The ester oxygen

        # Check each carbon adjacent to the carbonyl carbon if it begins a 10-carbon chain
        for neighbor in mol.GetAtomWithIdx(carbonyl_c_idx).GetNeighbors():
            if neighbor.GetAtomicNum() == 6: # carbon
                if neighbor.GetIdx() != oxygen_idx and find_carbon_chain(neighbor.GetIdx()):
                    return True, "Contains ester group with a decanoic acid (10-carbon) chain"

    return False, "Ester group without a proper decanoic acid chain"