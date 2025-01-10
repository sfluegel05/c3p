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
    
    # Define a SMARTS pattern for myo-inositol
    inositol_pattern = Chem.MolFromSmarts("[C@H]1(O)[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O")
    if not mol.HasSubstructMatch(inositol_pattern):
        return False, "No myo-inositol core structure found"
    
    # Define a SMARTS pattern for a phosphate group
    phosphate_pattern = Chem.MolFromSmarts("OP(O)(O)=O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) < 1:
        return False, "No phosphate groups found"
    
    return True, "Contains myo-inositol core with one or more phosphate groups"