"""
Classifies: CHEBI:25448 myo-inositol phosphate
"""
from rdkit import Chem

def is_myo_inositol_phosphate(smiles: str):
    """
    Determines if a molecule is a myo-inositol phosphate based on its SMILES string.
    A myo-inositol phosphate is an inositol phosphate in which the inositol component
    has myo-configuration.

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

    # Check for cyclohexane ring with myo-inositol stereochemistry
    myo_inositol_pattern = Chem.MolFromSmarts("[C@H]1([O])[C@@H]([O])[C@H]([O])[C@@H]([O])[C@H]([O])[C@@H]1[O]")
    if not mol.HasSubstructMatch(myo_inositol_pattern):
        return False, "No myo-inositol core structure found"
    
    # Check for presence of phosphate groups
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if len(phosphate_matches) == 0:
        return False, "No phosphate groups found"

    return True, "Contains myo-inositol core structure with phosphate groups"

# Testing the function
smiles_example = "O[C@H]1[C@H](OP(O)(O)=O)[C@H](O)[C@H](OP(O)(O)=O)[C@@H](O)[C@@H]1OP(O)(O)=O"
print(is_myo_inositol_phosphate(smiles_example))