"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    The molecule must contain a myo-inositol core with one or more phosphate groups directly attached.

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

    # Define a SMARTS pattern for the myo-inositol core, considering stereochemistry
    inositol_pattern = Chem.MolFromSmarts("O[C@@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")
    
    # Check if the molecule contains the myo-inositol core
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol core structure found"

    # Define a SMARTS pattern for phosphate groups, considering potential protonation states
    phosphate_pattern = Chem.MolFromSmarts("OP([O-])(=O)[O-]")
    
    # Check for multiple phosphate groups attached to the inositol core
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found directly attached to inositol core"

    # Check for large lipid-like or other complex molecules that are false positives
    # For instance, molecules with excessive carbons outside the inositol structure
    large_alkyl_chain_pattern = Chem.MolFromSmarts("C([C@@H](C)O)(C)C")
    if mol.HasSubstructMatch(large_alkyl_chain_pattern):
        return False, "False positive: large lipid or complex molecule detected"
    
    return True, "Contains a myo-inositol core with one or more phosphate groups attached directly"