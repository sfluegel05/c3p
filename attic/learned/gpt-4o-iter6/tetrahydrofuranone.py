"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is an oxolane with an oxo- substituent.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tetrahydrofuranone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Identifying a five-membered ring with one oxygen atom (oxolane)
    oxolane_pattern = Chem.MolFromSmarts("C1COCC1")  # 5-membered ring with O and 4 carbons
    if not mol.HasSubstructMatch(oxolane_pattern):
        return False, "No oxolane ring found"

    # Identifying the presence of a carbonyl group on the oxolane
    carbonyl_pattern = Chem.MolFromSmarts("C=O")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if not carbonyl_matches:
        return False, "No carbonyl group found on the ring"

    # Verify the carbonyl is attached to the ring
    for match in carbonyl_matches:
        # Check if any of the carbonyl carbon is part of the oxolane ring
        carbon = mol.GetAtomWithIdx(match[0])
        for neighbor in carbon.GetNeighbors():
            if neighbor.IsInRing() and neighbor.GetAtomicNum() == 8:  # Check if neighbor is oxygen in a ring
                return True, "Tetrahydrofuranone structure confirmed"

    return False, "Carbonyl group not part of an oxolane ring"