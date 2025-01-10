"""
Classifies: CHEBI:47016 tetrahydrofuranone
"""
from rdkit import Chem

def is_tetrahydrofuranone(smiles: str):
    """
    Determines if a molecule is a tetrahydrofuranone based on its SMILES string.
    A tetrahydrofuranone is any oxolane (tetrahydrofuran) having an oxo-substituent
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

    # Correct pattern for tetrahydrofuran (oxolane) ring [O]1CCCC1
    oxolane_pattern = Chem.MolFromSmarts("O1CCCC1")
    if not mol.HasSubstructMatch(oxolane_pattern):
        return False, "No oxolane (tetrahydrofuran) ring found"
    
    # Search for oxo group (C=O) positionally relaxed from ring but connected;
    oxo_relaxed = Chem.MolFromSmarts("[$(OC=O),$(OCC=O),$(CCC=O),$(CCCC=O)]")
    if not mol.HasSubstructMatch(oxo_relaxed):
        return False, "No oxo-substituent found on the tetrahydrofuran ring"

    return True, "Contains oxolane (tetrahydrofuran) ring with an oxo-substituent"