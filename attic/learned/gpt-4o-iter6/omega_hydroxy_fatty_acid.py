"""
Classifies: CHEBI:10615 omega-hydroxy fatty acid
"""
from rdkit import Chem

def is_omega_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an omega-hydroxy fatty acid based on its SMILES string.
    An omega-hydroxy fatty acid is a straight-chain fatty acid with:
    - A carboxyl group at position 1
    - A hydroxyl group at the omega (last) position

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an omega-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxyl group (-C(=O)O) at the end of the chain
    carboxyl_pattern = Chem.MolFromSmarts("C(=O)[OH]")
    if not mol.HasSubstructMatch(carboxyl_pattern):
        return False, "No terminal carboxyl group found"

    # Check for hydroxyl group (OH) at the omega position
    # A hydroxyl at the omega position would be at the end of the chain opposite to the carboxyl
    hydroxyl_pattern = Chem.MolFromSmarts("O")
    hydroxyl_matches = mol.GetSubstructMatches(hydroxyl_pattern)

    # Identify terminal hydroxyl which should mean at least two carbon separation from the carboxyl
    terminal_hydroxyl = False
    for match in hydroxyl_matches:
        atom = mol.GetAtomWithIdx(match[0])
        neighbors = [n for n in atom.GetNeighbors() if n.GetAtomicNum() == 6]
        if len(neighbors) == 1 and atom.GetDegree() == 1:  # Terminal OH check
            terminal_hydroxyl = True
            
    if not terminal_hydroxyl:
        return False, "No terminal omega-hydroxyl group found"

    # Ensure it is a straight-chain fatty acid (typically 7 or more carbons)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 7:
        return False, f"Chain length too short for fatty acid (found {carbon_count} carbons)"

    return True, "Matches omega-hydroxy fatty acid structure"