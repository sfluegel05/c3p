"""
Classifies: CHEBI:73754 thiosugar
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_thiosugar(smiles: str):
    """
    Determines if a molecule is a thiosugar based on its SMILES string.
    A thiosugar has a sugar-like cyclic structure where one of the oxygen atoms
    is replaced by a sulfur atom, often maintaining the hydroxyl configuration of sugars.

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
    
    # Define potential thiosugar patterns
    thiosugar_patterns = [
        Chem.MolFromSmarts("C1(O)[C@@H](O)[C@H](O)[C@H](O)[C@@H]1[S]"),  # General thiosugar pyranose pattern with sulfur
        Chem.MolFromSmarts("[C@H]1(O[C@H](S)[C@@H](O)[C@H](O)[C@H]1O)"),  # Thiosugar pattern with sulfur replacing oxygen
        Chem.MolFromSmarts("O[C@H]1[C@H](O)[C@@H](S)[C@H](O)[C@H](O1)")   # Another variant pattern
    ]

    # Check for thiosugar pattern match
    has_thiosugar = any(mol.HasSubstructMatch(pattern) for pattern in thiosugar_patterns)
    if not has_thiosugar:
        return False, "No thiosugar pattern found in the molecule"

    # Calculate potential sugar-like characteristics
    num_oh_groups = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'O' and any(n.GetSymbol() == 'H' for n in atom.GetNeighbors()))
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()

    if num_rings == 0:
        return False, "No ring structure found, unlikely to be a sugar"

    # Typically thiosugars are monosaccharides with one sulfur and 3-4 OH groups
    if num_oh_groups < 3:
        return False, f"Insufficient OH groups for a sugar-like structure, found only {num_oh_groups}"

    return True, "Molecule is a thiosugar: Shows thiosugar patterns and OH configuration"

__metadata__ = {
    'chemical_class': {   
        'name': 'thiosugar',
        'definition': 'Sugars in which one of the oxygen atoms is replaced by sulfur'
    }
}