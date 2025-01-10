"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is defined as any oxolane having an oxo-substituent
    at any position on the tetrahydrofuran ring.

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

    # Look for oxolane (tetrahydrofuran) ring pattern, [O]1CCCC1
    oxolane_pattern = Chem.MolFromSmarts("O1CCCCC1")
    if not mol.HasSubstructMatch(oxolane_pattern):
        return False, "No oxolane (tetrahydrofuran) ring found"
    
    # Look for oxo group (C=O) connected to the oxolane ring
    oxo_group_pattern = Chem.MolFromSmarts("C(=O)C1CCCCO1")
    if not mol.HasSubstructMatch(oxo_group_pattern):
        return False, "No oxo group attached to the oxolane ring found"

    return True, "Contains oxolane (tetrahydrofuran) ring with an oxo-substituent"