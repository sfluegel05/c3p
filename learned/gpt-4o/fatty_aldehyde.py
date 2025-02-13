"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is characterized by having an aldehyde group at one terminal of the carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty aldehyde, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string into a molecule object
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for the aldehyde carbonyl group pattern: R-C=O
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)[#6]")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No terminal aldehyde group (R-CHO) found"

    # Check each terminal aldehyde carbon to confirm it is part of a longer, acyclic chain
    for match in mol.GetSubstructMatches(aldehyde_pattern):
        aldehyde_carbon = mol.GetAtomWithIdx(match[0])

        # Ensure connected structure is an acyclic chain ideally with significant length
        neighbors = aldehyde_carbon.GetNeighbors()
        if len(neighbors) != 1:
            continue

        chain_start = neighbors[0]  # The C atom in R-CHO group that connects to the chain
        visited = set()
        carbon_count = 0
        to_visit = [chain_start]
        while to_visit:
            current = to_visit.pop()
            if current.GetIdx() in visited:
                continue
            visited.add(current.GetIdx())

            if current.GetAtomicNum() == 6:  # Count only carbon atoms
                carbon_count += 1
                for neighbor in current.GetNeighbors():
                    if neighbor.GetIdx() not in visited and not neighbor.IsInRing():
                        to_visit.append(neighbor)
        
        if carbon_count >= 5:
            return True, f"Has terminal aldehyde group in acyclic carbon chain of {carbon_count} carbons"

    return False, "Aldehyde group not terminal in proper chain structure"