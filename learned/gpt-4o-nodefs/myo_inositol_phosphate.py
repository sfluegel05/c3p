"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Broaden the SMARTS pattern for myo-inositol to include stereochemistry and varying hydroxyl positions
    inositol_patterns = [
        Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"),
        Chem.MolFromSmarts("[C@@H]1(O)[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)C1O")
    ]
    
    # Check if any of the myo-inositol patterns match
    inositol_match = any(mol.HasSubstructMatch(pattern) for pattern in inositol_patterns)
    if not inositol_match:
        return False, "No myo-inositol core structure found"
    
    # Define a SMARTS pattern for a phosphate group, considering variabilities in protonation
    phosphate_pattern = Chem.MolFromSmarts("O[P](=O)([O-])O")
    
    # Check for phosphate groups; must be directly attached to the inositol ring
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found directly attached to inositol core"

    # Further differentiation can be added here to avoid known false positive structures
    # For example, by checking that the molecule isn't a larger structure known to cause false positives
    
    return True, "Contains myo-inositol core with one or more phosphate groups attached directly"