"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar typically has a cyclic sugar structure where one or more oxygen atoms
    are replaced by sulfur atoms, maintaining a sugar-like hydroxyl pattern.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a thiosugar, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define thiosugar pattern where sulfur replaces an oxygen in the sugar cycle
    patterns = [
        Chem.MolFromSmarts("[C@H]1([C@H](O)[C@@H](S)[C@H](O)[C@H]1O)"),  # Thiosugar pattern with sulfur in the ring
        Chem.MolFromSmarts("[C@@H]1([O][C@H](S)[C@H](O)[C@H](O)[C@H]1O)"),  # Variant pattern
    ]

    # Check for thiosugar pattern match
    has_thiosugar_structure = any(mol.HasSubstructMatch(pattern) for pattern in patterns)
    
    if not has_thiosugar_structure:
        return False, "No thiosugar pattern found in the molecule"
    
    # Calculate sugar-like characteristics
    num_oh_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(n.GetSymbol() == 'H' for n in atom.GetNeighbors()))
    num_sulfur_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'S')
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    # Verify the presence of ring and OH groups for sugar-like structure
    if num_rings == 0:
        return False, "No ring structure found, unlikely to be a sugar"
    if num_oh_groups < 3:
        return False, f"Insufficient OH groups for a sugar-like structure, found only {num_oh_groups}"
    if num_sulfur_atoms < 1:
        return False, "No sulfur found replacing an oxygen atom in the structure"

    return True, "Molecule is a thiosugar: Shows thiosugar patterns and OH configuration"

__metadata__ = {
    'chemical_class': {   
        'name': 'thiosugar',
        'definition': 'Sugars where one or more oxygen atoms are replaced by sulfur'
    }
}