"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate with a specific stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a myo-inositol phosphate, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS pattern for the cyclohexane ring typical in myo-inositol
    # It should have all OH groups and flexible stereochemistry captured
    myo_inositol_pattern = Chem.MolFromSmarts("C1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O")

    # Check for the myo-inositol backbone with possible stereochemical variations
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol backbone or mismatched stereochemistry found"
    
    # Define SMARTS pattern for phosphate groups, including different possible forms
    phosphate_pattern = Chem.MolFromSmarts("OP(=O)(O)O")

    # Check for presence of phosphate groups
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"
    
    return True, "Contains myo-inositol backbone with one or more phosphate groups"