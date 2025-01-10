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

    # Correct pattern for a tetrahydrofuran (oxolane) ring - [O]1CCCC1
    oxolane_pattern = Chem.MolFromSmarts("O1CCCC1")
    if not mol.HasSubstructMatch(oxolane_pattern):
        return False, "No oxolane (tetrahydrofuran) ring found"
    
    # Search for oxo group (C=O) attached directly to the oxolane ring
    # Modification: Looking at carbonyls bonded to any carbon on the ring
    # more predictively captures the expected structures
    oxo_attached_to_ring = Chem.MolFromSmarts("O1CCCC1C=O")
    if not mol.HasSubstructMatch(oxo_attached_to_ring):
        # As a fallback, expand the view to check for nearby
        nearby_oxo_pattern = Chem.MolFromSmarts("C=O")
        attached_atoms = mol.GetSubstructMatch(oxolane_pattern)
        for atom_id in attached_atoms:
            atom = mol.GetAtomWithIdx(atom_id)
            for neighbor in atom.GetNeighbors():
                if neighbor.HasQueryMatch(nearby_oxo_pattern):
                    return True, "Contains oxolane (tetrahydrofuran) ring with nearby oxo-substituent"
        return False, "No oxo-substituent found on the tetrahydrofuran ring appropriately connected"
    
    return True, "Contains oxolane (tetrahydrofuran) ring with an oxo-substituent"