"""
Classifies: CHEBI:35746 fatty aldehyde
"""
from rdkit import Chem

def is_fatty_aldehyde(smiles: str):
    """
    Determines if a molecule is a fatty aldehyde based on its SMILES string.
    A fatty aldehyde is characterized by having an aldehyde group at the end of a carbon chain.

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

    # Look for the aldehyde carbonyl group pattern: C=O
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)C")
    if not mol.HasSubstructMatch(aldehyde_pattern):
        return False, "No aldehyde group found"

    # Check each terminal aldehyde carbon to confirm it is part of a longer, acyclic chain
    for match in mol.GetSubstructMatches(aldehyde_pattern):
        carbon_atom = mol.GetAtomWithIdx(match[0])
        if carbon_atom.GetTotalDegree() == 2:
            # Ensure the connected chain is acyclic and reasonably long (4+ carbons)
            # Walk through the chain and count consecutive carbon bonds
            visited = set()
            carbon_count = 0
            to_visit = [carbon_atom]
            while to_visit:
                current = to_visit.pop()
                visited.add(current.GetIdx())
                if current.GetAtomicNum() == 6:  # Count only carbon atoms
                    carbon_count += 1
                    for neighbor in current.GetNeighbors():
                        if neighbor.GetIdx() not in visited and not neighbor.IsInRing():
                            to_visit.append(neighbor)
            if carbon_count >= 4:
                return True, f"Has terminal aldehyde group in acyclic carbon chain of {carbon_count} carbons"

    return False, "Aldehyde group not terminal in proper chain structure"